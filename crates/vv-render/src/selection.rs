//! Selection outline: a screen-space halo around every visible pixel
//! whose pick id is selected (`shaders/selection_outline.wgsl`). It reads
//! the frame's own id target, so it follows whatever representation drew
//! the atom, is occluded exactly as the atom is, and moves with it in the
//! same frame. Drawn onto the display only, never into `read_color`.

use bytemuck::{Pod, Zeroable};
use wgpu::util::DeviceExt;

use crate::renderer::{storage_entry, texture_entry, uniform_entry, PickRange, COLOR_FORMAT};

/// The selected pick ids of one pick range, indexed like `Pick::Atom {
/// item }`: `render_all`'s items, then cartoons, then patch surfaces.
/// Bit `i` of word `i / 32` (LSB first) marks local id `i`: an atom (or
/// cartoon section, or tube sample) in `atoms`, a bond in `bonds`. Bits
/// past the range's own count are ignored; missing words are unselected.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ItemSelection {
    pub atoms: Vec<u32>,
    pub bonds: Vec<u32>,
}

impl ItemSelection {
    /// `n` bits, bit `i` = `selected(i)`.
    pub fn bits(n: usize, selected: impl Fn(usize) -> bool) -> Vec<u32> {
        let mut words = vec![0u32; n.div_ceil(32)];
        for i in (0..n).filter(|&i| selected(i)) {
            words[i / 32] |= 1 << (i % 32);
        }
        words
    }
}

/// Layout must match `Style` in `shaders/selection_outline.wgsl`.
#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct StyleUniform {
    color: [f32; 4],
    width: u32,
    _pad: [u32; 3],
}

pub(crate) struct SelectionOutline {
    pipeline: wgpu::RenderPipeline,
    ids_layout: wgpu::BindGroupLayout,
    bits_layout: wgpu::BindGroupLayout,
    style: wgpu::Buffer,
    atom_bits: wgpu::Buffer,
    bond_bits: wgpu::Buffer,
    ids_group: wgpu::BindGroup,
    bits_group: wgpu::BindGroup,
    items: Vec<ItemSelection>,
    /// The ranges the uploaded bits were packed for; `None` forces a pack.
    packed_for: Option<Vec<PickRange>>,
    /// Whether any uploaded bit is set; the pass is skipped otherwise.
    any: bool,
}

impl SelectionOutline {
    pub fn new(device: &wgpu::Device, fullscreen: &str, id_view: &wgpu::TextureView) -> Self {
        let fragment = wgpu::ShaderStages::FRAGMENT;
        let ids_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("selection ids"),
            entries: &[texture_entry(0, wgpu::TextureSampleType::Uint, fragment)],
        });
        let bits_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("selection bits"),
            entries: &[
                storage_entry(0, true, fragment),
                storage_entry(1, true, fragment),
                uniform_entry(2, fragment),
            ],
        });
        let style = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("selection style"),
            size: std::mem::size_of::<StyleUniform>() as u64,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        let atom_bits = bits_buffer(device, "selection atom bits", &[0]);
        let bond_bits = bits_buffer(device, "selection bond bits", &[0]);
        let bits_group = bits_group(device, &bits_layout, &atom_bits, &bond_bits, &style);
        let ids_group = ids_group(device, &ids_layout, id_view);
        let pipeline = pipeline(device, fullscreen, &ids_layout, &bits_layout);
        Self {
            pipeline,
            ids_layout,
            bits_layout,
            style,
            atom_bits,
            bond_bits,
            ids_group,
            bits_group,
            items: Vec::new(),
            packed_for: None,
            any: false,
        }
    }

    pub fn set(&mut self, items: Vec<ItemSelection>) {
        if items != self.items {
            self.items = items;
            self.packed_for = None;
        }
    }

    /// Call when the id target is recreated.
    pub fn rebind_ids(&mut self, device: &wgpu::Device, id_view: &wgpu::TextureView) {
        self.ids_group = ids_group(device, &self.ids_layout, id_view);
    }

    /// Re-packs and uploads the bits only when the selection or the id
    /// ranges changed since the last upload.
    pub fn prepare(
        &mut self,
        device: &wgpu::Device,
        queue: &wgpu::Queue,
        ranges: &[PickRange],
        color: [f32; 4],
        width: u32,
    ) {
        let style = StyleUniform {
            color,
            width: width.max(1),
            _pad: [0; 3],
        };
        queue.write_buffer(&self.style, 0, bytemuck::bytes_of(&style));
        if self.packed_for.as_deref() == Some(ranges) {
            return;
        }
        let (atoms, bonds) = pack(&self.items, ranges);
        self.any = atoms.iter().chain(&bonds).any(|&w| w != 0);
        let atoms_grown = self.upload(device, queue, &atoms, true);
        let bonds_grown = self.upload(device, queue, &bonds, false);
        if atoms_grown || bonds_grown {
            self.bits_group = bits_group(
                device,
                &self.bits_layout,
                &self.atom_bits,
                &self.bond_bits,
                &self.style,
            );
        }
        self.packed_for = Some(ranges.to_vec());
    }

    /// Writes `words` over the whole buffer (zero-padded, so no stale bits
    /// survive past them), reallocating when it outgrew it.
    fn upload(
        &mut self,
        device: &wgpu::Device,
        queue: &wgpu::Queue,
        words: &[u32],
        atoms: bool,
    ) -> bool {
        let buffer = if atoms {
            &mut self.atom_bits
        } else {
            &mut self.bond_bits
        };
        let capacity = (buffer.size() / 4) as usize;
        if words.len() > capacity {
            *buffer = bits_buffer(device, buffer_label(atoms), words);
            return true;
        }
        let mut padded = words.to_vec();
        padded.resize(capacity, 0);
        queue.write_buffer(buffer, 0, bytemuck::cast_slice(&padded));
        false
    }

    /// Draws the outline over `target` (the display image) in place.
    pub fn draw(&self, encoder: &mut wgpu::CommandEncoder, target: &wgpu::TextureView) {
        if !self.any {
            return;
        }
        let mut pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
            label: Some("selection outline"),
            color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                view: target,
                resolve_target: None,
                ops: wgpu::Operations {
                    load: wgpu::LoadOp::Load,
                    store: wgpu::StoreOp::Store,
                },
                depth_slice: None,
            })],
            ..Default::default()
        });
        pass.set_pipeline(&self.pipeline);
        pass.set_bind_group(0, &self.ids_group, &[]);
        pass.set_bind_group(1, &self.bits_group, &[]);
        pass.draw(0..3, 0..1);
    }
}

fn buffer_label(atoms: bool) -> &'static str {
    if atoms {
        "selection atom bits"
    } else {
        "selection bond bits"
    }
}

/// Frame-wide atom and bond bitsets: each item's local bits shifted to
/// its range's base, clipped to its count.
fn pack(items: &[ItemSelection], ranges: &[PickRange]) -> (Vec<u32>, Vec<u32>) {
    let atom_total = ranges.iter().map(|r| r.atom_base + r.atom_count).max();
    let bond_total = ranges.iter().map(|r| r.bond_base + r.bond_count).max();
    let mut atoms = vec![0u32; atom_total.unwrap_or(0).div_ceil(32) as usize];
    let mut bonds = vec![0u32; bond_total.unwrap_or(0).div_ceil(32) as usize];
    for (item, range) in items.iter().zip(ranges) {
        place(&mut atoms, &item.atoms, range.atom_base, range.atom_count);
        place(&mut bonds, &item.bonds, range.bond_base, range.bond_count);
    }
    (atoms, bonds)
}

fn place(out: &mut [u32], local: &[u32], base: u32, count: u32) {
    for (w, &word) in local.iter().enumerate() {
        let mut word = word;
        while word != 0 {
            let i = w as u32 * 32 + word.trailing_zeros();
            word &= word - 1;
            if i >= count {
                return;
            }
            let g = base + i;
            out[(g / 32) as usize] |= 1 << (g % 32);
        }
    }
}

fn bits_buffer(device: &wgpu::Device, label: &str, words: &[u32]) -> wgpu::Buffer {
    device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some(label),
        contents: bytemuck::cast_slice(words),
        usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
    })
}

fn ids_group(
    device: &wgpu::Device,
    layout: &wgpu::BindGroupLayout,
    id_view: &wgpu::TextureView,
) -> wgpu::BindGroup {
    device.create_bind_group(&wgpu::BindGroupDescriptor {
        label: Some("selection ids"),
        layout,
        entries: &[wgpu::BindGroupEntry {
            binding: 0,
            resource: wgpu::BindingResource::TextureView(id_view),
        }],
    })
}

fn bits_group(
    device: &wgpu::Device,
    layout: &wgpu::BindGroupLayout,
    atom_bits: &wgpu::Buffer,
    bond_bits: &wgpu::Buffer,
    style: &wgpu::Buffer,
) -> wgpu::BindGroup {
    device.create_bind_group(&wgpu::BindGroupDescriptor {
        label: Some("selection bits"),
        layout,
        entries: &[
            wgpu::BindGroupEntry {
                binding: 0,
                resource: atom_bits.as_entire_binding(),
            },
            wgpu::BindGroupEntry {
                binding: 1,
                resource: bond_bits.as_entire_binding(),
            },
            wgpu::BindGroupEntry {
                binding: 2,
                resource: style.as_entire_binding(),
            },
        ],
    })
}

fn pipeline(
    device: &wgpu::Device,
    fullscreen: &str,
    ids_layout: &wgpu::BindGroupLayout,
    bits_layout: &wgpu::BindGroupLayout,
) -> wgpu::RenderPipeline {
    let module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
        label: Some("selection outline"),
        source: wgpu::ShaderSource::Wgsl(
            format!(
                "{fullscreen}\n{}",
                include_str!("shaders/selection_outline.wgsl")
            )
            .into(),
        ),
    });
    let layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
        label: Some("selection outline"),
        bind_group_layouts: &[Some(ids_layout), Some(bits_layout)],
        ..Default::default()
    });
    device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
        label: Some("selection outline"),
        layout: Some(&layout),
        vertex: wgpu::VertexState {
            module: &module,
            entry_point: Some("vs_fullscreen"),
            compilation_options: Default::default(),
            buffers: &[],
        },
        primitive: Default::default(),
        depth_stencil: None,
        multisample: Default::default(),
        fragment: Some(wgpu::FragmentState {
            module: &module,
            entry_point: Some("fs_selection_outline"),
            compilation_options: Default::default(),
            targets: &[Some(wgpu::ColorTargetState {
                format: COLOR_FORMAT,
                blend: None,
                write_mask: wgpu::ColorWrites::ALL,
            })],
        }),
        multiview_mask: None,
        cache: None,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pack_shifts_each_item_to_its_base_and_clips_to_its_count() {
        let ranges = [
            PickRange {
                atom_base: 0,
                atom_count: 3,
                bond_base: 0,
                bond_count: 2,
            },
            PickRange {
                atom_base: 3,
                atom_count: 40,
                bond_base: 2,
                bond_count: 1,
            },
        ];
        let items = [
            ItemSelection {
                atoms: ItemSelection::bits(5, |i| i == 1 || i == 4),
                bonds: vec![],
            },
            ItemSelection {
                atoms: ItemSelection::bits(40, |i| i == 0 || i == 33),
                bonds: ItemSelection::bits(1, |_| true),
            },
        ];
        let (atoms, bonds) = pack(&items, &ranges);
        let set: Vec<u32> = (0..atoms.len() as u32 * 32)
            .filter(|&g| atoms[(g / 32) as usize] & (1 << (g % 32)) != 0)
            .collect();
        assert_eq!(set, [1, 3, 36]);
        assert_eq!(bonds, [0b100]);
    }
}
