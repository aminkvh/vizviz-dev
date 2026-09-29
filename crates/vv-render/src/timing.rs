//! GPU timestamp queries with a small ring of readback buffers so reading
//! results never stalls the frame that produced them.
//!
//! Each slot moves free -> copied -> mapping -> ready -> free. A slot whose
//! mapping has not completed yet (the GPU is several frames behind) is
//! skipped for that frame rather than written into while mapped.

use std::sync::atomic::{AtomicU8, Ordering};
use std::sync::Arc;

use crate::context::GpuContext;

const SLOTS: usize = 3;
/// Timestamps per frame: cull, draw, depth-pyramid and post begin/end
/// pairs. The first and last also bound the whole frame.
pub const QUERIES: u32 = 8;

const FREE: u8 = 0;
const COPIED: u8 = 1;
const MAPPING: u8 = 2;
const READY: u8 = 3;

pub struct GpuTimer {
    query_set: wgpu::QuerySet,
    resolve: wgpu::Buffer,
    readback: Vec<(wgpu::Buffer, Arc<AtomicU8>)>,
    period_ns: f32,
    frame: usize,
    /// Slot whose copy was encoded this frame and must be mapped after submit.
    submitted: Option<usize>,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct PassTimes {
    /// Phase-1 cull.
    pub cull_ms: f32,
    /// Phase-1 opaque pass.
    pub draw_ms: f32,
    /// The final depth pyramid.
    pub hiz_ms: f32,
    /// AO, composite and FXAA.
    pub post_ms: f32,
    /// First pass start to last pass end: everything above plus the
    /// occlusion phase 2 (mid-frame pyramid, cull, draw).
    pub frame_ms: f32,
}

impl PassTimes {
    pub fn total_ms(&self) -> f32 {
        self.frame_ms
    }
}

impl GpuTimer {
    pub fn new(ctx: &GpuContext) -> Option<Self> {
        if !ctx.timestamps {
            return None;
        }
        let device = &ctx.device;
        let query_set = device.create_query_set(&wgpu::QuerySetDescriptor {
            label: Some("frame timestamps"),
            ty: wgpu::QueryType::Timestamp,
            count: QUERIES,
        });
        let size = QUERIES as u64 * 8;
        let resolve = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("timestamp resolve"),
            size,
            usage: wgpu::BufferUsages::QUERY_RESOLVE | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: false,
        });
        let readback = (0..SLOTS)
            .map(|i| {
                let buffer = device.create_buffer(&wgpu::BufferDescriptor {
                    label: Some(&format!("timestamp readback {i}")),
                    size,
                    usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
                    mapped_at_creation: false,
                });
                (buffer, Arc::new(AtomicU8::new(FREE)))
            })
            .collect();
        Some(Self {
            query_set,
            resolve,
            readback,
            period_ns: ctx.queue.get_timestamp_period(),
            frame: 0,
            submitted: None,
        })
    }

    pub fn query_set(&self) -> &wgpu::QuerySet {
        &self.query_set
    }

    /// Call once per frame after all passes were recorded. Harvests the
    /// oldest completed frame's times (if any) and, when this frame's slot
    /// is free, schedules this frame's queries for readback.
    pub fn end_frame(&mut self, encoder: &mut wgpu::CommandEncoder) -> Option<PassTimes> {
        let slot = self.frame % SLOTS;
        self.frame += 1;
        let (buffer, state) = &self.readback[slot];

        let mut result = None;
        if state.load(Ordering::Acquire) == READY {
            let data = buffer.slice(..).get_mapped_range().expect("mapped range");
            let stamps: &[u64] = bytemuck::cast_slice(&data);
            let ms = |a: u64, b: u64| b.saturating_sub(a) as f32 * self.period_ns * 1e-6;
            result = Some(PassTimes {
                cull_ms: ms(stamps[0], stamps[1]),
                draw_ms: ms(stamps[2], stamps[3]),
                hiz_ms: ms(stamps[4], stamps[5]),
                post_ms: ms(stamps[6], stamps[7]),
                frame_ms: ms(stamps[0], stamps[7]),
            });
            drop(data);
            buffer.unmap();
            state.store(FREE, Ordering::Release);
        }

        if state.load(Ordering::Acquire) == FREE {
            encoder.resolve_query_set(&self.query_set, 0..QUERIES, &self.resolve, 0);
            encoder.copy_buffer_to_buffer(&self.resolve, 0, buffer, 0, QUERIES as u64 * 8);
            state.store(COPIED, Ordering::Release);
            self.submitted = Some(slot);
        }
        result
    }

    /// Call right after the frame's command buffer is submitted.
    pub fn after_submit(&mut self, ctx: &GpuContext) {
        if let Some(slot) = self.submitted.take() {
            let (buffer, state) = &self.readback[slot];
            state.store(MAPPING, Ordering::Release);
            let flag = state.clone();
            buffer.slice(..).map_async(wgpu::MapMode::Read, move |r| {
                flag.store(if r.is_ok() { READY } else { FREE }, Ordering::Release);
            });
        }
        // Drive completed map callbacks without blocking.
        let _ = ctx.device.poll(wgpu::PollType::Poll);
    }
}
