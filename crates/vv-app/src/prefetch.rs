//! Reads a streamed trajectory's next frames on worker threads while it
//! plays, so the frame the timeline steps to is already in memory: the
//! UI thread never waits on a decode (an XTC frame of a million atoms is
//! tens of milliseconds). Frames land in the structure's own cache
//! (`vv_core::Structure::try_frame`).

use std::collections::HashSet;
use std::sync::{Arc, Mutex};

use vv_core::Structure;
use vv_scene::StructureId;

use crate::gpu_cache::Waker;

#[derive(Default)]
pub struct Prefetch {
    in_flight: Arc<Mutex<HashSet<(StructureId, usize)>>>,
}

impl Prefetch {
    /// Starts reading the `ahead` frames from `from` on (wrapping round
    /// when `looping`) that are neither in memory nor already coming;
    /// `wake` runs as each arrives.
    pub fn ahead(
        &self,
        id: StructureId,
        structure: &Structure,
        from: usize,
        ahead: usize,
        looping: bool,
        wake: Option<Waker>,
    ) {
        if !structure.is_streamed() {
            return;
        }
        let count = structure.frame_count();
        for step in 0..ahead.min(count) {
            let frame = from + step;
            let frame = match (frame < count, looping) {
                (true, _) => frame,
                (false, true) => frame % count,
                (false, false) => break,
            };
            if structure.has_frame(frame)
                || !self
                    .in_flight
                    .lock()
                    .expect("prefetch set poisoned")
                    .insert((id, frame))
            {
                continue;
            }
            let (structure, in_flight, wake) =
                (structure.clone(), self.in_flight.clone(), wake.clone());
            rayon::spawn(move || {
                // An unreadable frame is reported where it is shown.
                let _ = structure.try_frame(frame);
                in_flight
                    .lock()
                    .expect("prefetch set poisoned")
                    .remove(&(id, frame));
                if let Some(wake) = wake {
                    wake();
                }
            });
        }
    }
}
