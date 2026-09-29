//! A structure = one topology + one or more coordinate sets.
//!
//! Both halves are immutable behind `Arc`: edits build a new value and swap
//! the pointer, so GPU uploads, undo history, and Python views can keep
//! holding the old one safely.
//!
//! Frames are either all in memory (a parsed file's models) or streamed
//! from a trajectory file ([`FrameSource`]) through a cache capped in
//! bytes, so a trajectory larger than memory plays: frames are read when
//! asked for, and the least recently used ones are dropped.

use std::collections::HashMap;
use std::fmt;
use std::sync::{Arc, Mutex};

use glam::Vec3;

use crate::{CoordSet, Topology};

/// Where a streamed structure's frames come from: a trajectory file read
/// on demand. Called from any thread.
pub trait FrameSource: Send + Sync {
    fn frame_count(&self) -> usize;
    /// Frame `index`'s positions, one per atom.
    fn read(&self, index: usize) -> Result<Vec<Vec3>, String>;
}

#[derive(Clone, Debug)]
pub struct Structure {
    pub topology: Arc<Topology>,
    frames: Frames,
}

#[derive(Clone)]
enum Frames {
    Loaded(Vec<Arc<CoordSet>>),
    Streamed(Arc<FrameStream>),
}

impl fmt::Debug for Frames {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Frames::Loaded(frames) => write!(f, "Loaded({} frames)", frames.len()),
            Frames::Streamed(s) => write!(f, "Streamed({} frames)", s.source.frame_count()),
        }
    }
}

struct FrameStream {
    source: Box<dyn FrameSource>,
    atom_count: usize,
    /// Frame 0 stays resident: it is what a failed read falls back to.
    first: Arc<CoordSet>,
    cache: Mutex<FrameCache>,
}

#[derive(Default)]
struct FrameCache {
    frames: HashMap<usize, (Arc<CoordSet>, u64)>,
    clock: u64,
    budget: usize,
    /// The last read that failed, for the app to report.
    error: Option<String>,
}

impl FrameCache {
    fn frame_bytes(atoms: usize) -> usize {
        atoms * std::mem::size_of::<Vec3>()
    }

    /// Adds `frame`, dropping the least recently used ones past the budget.
    fn insert(&mut self, index: usize, frame: Arc<CoordSet>) {
        self.clock += 1;
        let bytes = Self::frame_bytes(frame.len());
        let keep = (self.budget / bytes.max(1)).max(2);
        self.frames.insert(index, (frame, self.clock));
        while self.frames.len() > keep {
            let oldest = self
                .frames
                .iter()
                .min_by_key(|(_, (_, used))| *used)
                .map(|(&k, _)| k)
                .expect("over budget, so not empty");
            self.frames.remove(&oldest);
        }
    }
}

#[derive(thiserror::Error, Debug, PartialEq, Eq)]
pub enum StructureError {
    #[error("coordinate set {index} has {actual} positions, topology has {expected} atoms")]
    CoordCount {
        index: usize,
        expected: usize,
        actual: usize,
    },
    #[error("a structure needs at least one coordinate set")]
    NoCoords,
    #[error("reading frame {index}: {message}")]
    Read { index: usize, message: String },
}

/// Bytes of frames a streamed structure keeps in memory by default.
pub const DEFAULT_FRAME_BUDGET: usize = 1 << 30;

impl Structure {
    pub fn new(topology: Topology, coords: CoordSet) -> Result<Self, StructureError> {
        Self::with_frames(topology, vec![coords])
    }

    pub fn with_frames(
        mut topology: Topology,
        frames: Vec<CoordSet>,
    ) -> Result<Self, StructureError> {
        if frames.is_empty() {
            return Err(StructureError::NoCoords);
        }
        let expected = topology.atom_count();
        for (index, f) in frames.iter().enumerate() {
            check_count(index, f.len(), expected)?;
        }
        topology.assign_residue_classes(frames[0].positions());
        Ok(Self {
            topology: Arc::new(topology),
            frames: Frames::Loaded(frames.into_iter().map(Arc::new).collect()),
        })
    }

    /// Frames read from `source` on demand, keeping up to `budget` bytes
    /// of them (plus frame 0) in memory.
    pub fn streamed(
        mut topology: Arc<Topology>,
        source: Box<dyn FrameSource>,
        budget: usize,
    ) -> Result<Self, StructureError> {
        if source.frame_count() == 0 {
            return Err(StructureError::NoCoords);
        }
        let expected = topology.atom_count();
        let first = source
            .read(0)
            .map_err(|message| StructureError::Read { index: 0, message })?;
        check_count(0, first.len(), expected)?;
        if topology.residue_class.len() != topology.residue_count() {
            Arc::make_mut(&mut topology).assign_residue_classes(&first);
        }
        Ok(Self {
            topology,
            frames: Frames::Streamed(Arc::new(FrameStream {
                source,
                atom_count: expected,
                first: Arc::new(CoordSet::new(first)),
                cache: Mutex::new(FrameCache {
                    budget,
                    ..FrameCache::default()
                }),
            })),
        })
    }

    pub fn atom_count(&self) -> usize {
        self.topology.atom_count()
    }

    pub fn frame_count(&self) -> usize {
        match &self.frames {
            Frames::Loaded(frames) => frames.len(),
            Frames::Streamed(s) => s.source.frame_count(),
        }
    }

    /// Identifies this coordinate source: equal only for structures
    /// sharing the same frames, so a swapped-in trajectory is detectable.
    pub fn frames_id(&self) -> usize {
        match &self.frames {
            Frames::Loaded(frames) => Arc::as_ptr(&frames[0]) as usize,
            Frames::Streamed(s) => Arc::as_ptr(s) as usize,
        }
    }

    pub fn is_streamed(&self) -> bool {
        matches!(self.frames, Frames::Streamed(_))
    }

    /// Frame `index`, reading it if it is streamed and not in memory. For
    /// analysis: an unreadable frame is an error, never a substitute.
    pub fn try_frame(&self, index: usize) -> Result<Arc<CoordSet>, StructureError> {
        let s = match &self.frames {
            Frames::Loaded(frames) => return Ok(frames[index].clone()),
            Frames::Streamed(s) => s,
        };
        if index == 0 {
            return Ok(s.first.clone());
        }
        {
            let mut cache = s.cache.lock().expect("frame cache poisoned");
            cache.clock += 1;
            let clock = cache.clock;
            if let Some((frame, used)) = cache.frames.get_mut(&index) {
                *used = clock;
                return Ok(frame.clone());
            }
        }
        // Read without the lock, so other threads' frames aren't held up.
        let read = s
            .source
            .read(index)
            .map_err(|message| StructureError::Read { index, message })
            .and_then(|p| {
                check_count(index, p.len(), s.atom_count)?;
                Ok(Arc::new(CoordSet::new(p)))
            });
        let mut cache = s.cache.lock().expect("frame cache poisoned");
        match read {
            Ok(frame) => {
                cache.insert(index, frame.clone());
                Ok(frame)
            }
            Err(e) => {
                cache.error = Some(e.to_string());
                Err(e)
            }
        }
    }

    /// Frame `index`, for drawing: an unreadable streamed frame shows as
    /// frame 0 (see [`Structure::take_read_error`]).
    pub fn frame(&self, index: usize) -> Arc<CoordSet> {
        self.try_frame(index)
            .unwrap_or_else(|_| match &self.frames {
                Frames::Loaded(frames) => frames[0].clone(),
                Frames::Streamed(s) => s.first.clone(),
            })
    }

    /// Frame `index` by reference, when it lives as long as the structure:
    /// every in-memory frame, and a streamed structure's frame 0. `None`
    /// for other streamed frames, which may be dropped from memory.
    pub fn resident_frame(&self, index: usize) -> Option<&CoordSet> {
        match &self.frames {
            Frames::Loaded(frames) => frames.get(index).map(|f| &**f),
            Frames::Streamed(s) => (index == 0).then_some(&*s.first),
        }
    }

    /// Whether frame `index` is in memory now (reading it would not block).
    pub fn has_frame(&self, index: usize) -> bool {
        match &self.frames {
            Frames::Loaded(_) => true,
            Frames::Streamed(s) => {
                index == 0
                    || s.cache
                        .lock()
                        .expect("frame cache poisoned")
                        .frames
                        .contains_key(&index)
            }
        }
    }

    /// The last streamed frame that could not be read, once.
    pub fn take_read_error(&self) -> Option<String> {
        match &self.frames {
            Frames::Loaded(_) => None,
            Frames::Streamed(s) => s.cache.lock().expect("frame cache poisoned").error.take(),
        }
    }
}

fn check_count(index: usize, actual: usize, expected: usize) -> Result<(), StructureError> {
    if actual == expected {
        Ok(())
    } else {
        Err(StructureError::CoordCount {
            index,
            expected,
            actual,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::glam::Vec3;
    use crate::Element;

    #[test]
    fn rejects_mismatched_coordinate_count() {
        let topology = Topology {
            element: vec![Element::CARBON; 3],
            ..Default::default()
        };
        let err = Structure::new(topology, CoordSet::new(vec![Vec3::ZERO; 2])).unwrap_err();
        assert_eq!(
            err,
            StructureError::CoordCount {
                index: 0,
                expected: 3,
                actual: 2
            }
        );
    }

    #[test]
    fn rejects_no_frames() {
        assert_eq!(
            Structure::with_frames(Topology::default(), vec![]).unwrap_err(),
            StructureError::NoCoords
        );
    }

    /// Frame `k` has every atom at `(k, 0, 0)`; frame 3 can't be read.
    struct Counting {
        reads: std::sync::atomic::AtomicUsize,
    }

    impl FrameSource for Counting {
        fn frame_count(&self) -> usize {
            5
        }
        fn read(&self, index: usize) -> Result<Vec<Vec3>, String> {
            self.reads
                .fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            if index == 3 {
                return Err("bad frame".into());
            }
            Ok(vec![Vec3::new(index as f32, 0.0, 0.0); 2])
        }
    }

    fn streamed(budget: usize) -> (Structure, Arc<Counting>) {
        let topology = Arc::new(Topology {
            element: vec![Element::CARBON; 2],
            ..Default::default()
        });
        let source = Arc::new(Counting {
            reads: Default::default(),
        });
        struct Shared(Arc<Counting>);
        impl FrameSource for Shared {
            fn frame_count(&self) -> usize {
                self.0.frame_count()
            }
            fn read(&self, index: usize) -> Result<Vec<Vec3>, String> {
                self.0.read(index)
            }
        }
        let s = Structure::streamed(topology, Box::new(Shared(source.clone())), budget).unwrap();
        (s, source)
    }

    #[test]
    fn streamed_frames_are_read_once_while_they_fit() {
        let (s, source) = streamed(1 << 20);
        assert_eq!(s.frame_count(), 5);
        assert_eq!(s.frame(2).positions()[0].x, 2.0);
        assert_eq!(s.frame(2).positions()[0].x, 2.0);
        // Frame 0 at open, frame 2 once.
        assert_eq!(source.reads.load(std::sync::atomic::Ordering::Relaxed), 2);
        assert!(s.has_frame(2) && !s.has_frame(4));
    }

    #[test]
    fn the_least_recently_used_frame_goes_first_past_the_budget() {
        // Room for two frames of two atoms.
        let (s, _) = streamed(2 * 2 * std::mem::size_of::<Vec3>());
        s.frame(1);
        s.frame(2);
        s.frame(1);
        s.frame(4);
        assert!(s.has_frame(1) && s.has_frame(4) && !s.has_frame(2));
    }

    #[test]
    fn an_unreadable_frame_is_an_error_for_analysis_and_frame_zero_for_drawing() {
        let (s, _) = streamed(1 << 20);
        assert!(matches!(
            s.try_frame(3),
            Err(StructureError::Read { index: 3, .. })
        ));
        assert_eq!(s.frame(3).positions()[0].x, 0.0);
        assert!(s.take_read_error().is_some());
        assert!(s.take_read_error().is_none());
    }

    #[test]
    fn frames_share_one_topology() {
        let topology = Topology {
            element: vec![Element::CARBON; 2],
            ..Default::default()
        };
        let s = Structure::with_frames(
            topology,
            vec![
                CoordSet::new(vec![Vec3::ZERO; 2]),
                CoordSet::new(vec![Vec3::ONE; 2]),
            ],
        )
        .unwrap();
        assert_eq!(s.atom_count(), 2);
        assert_eq!(s.frame_count(), 2);
        assert_eq!(s.frame(1).positions()[0], Vec3::ONE);
    }
}
