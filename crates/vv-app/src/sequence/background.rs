//! Work that runs off the UI thread and is polled from the frame loop: a
//! result per key, computed once per stamp on a worker thread that wakes
//! the UI when it finishes.

use std::collections::HashMap;
use std::hash::Hash;
use std::sync::mpsc::{channel, Receiver, TryRecvError};
use std::sync::Arc;

enum State<T> {
    Running(Receiver<T>),
    Done(Arc<T>),
    /// The worker died without a result; not retried for this stamp.
    Failed,
}

struct Slot<S, T> {
    stamp: S,
    state: State<T>,
}

pub struct Jobs<K, S, T> {
    slots: HashMap<K, Slot<S, T>>,
}

impl<K, S, T> Default for Jobs<K, S, T> {
    fn default() -> Self {
        Self {
            slots: HashMap::new(),
        }
    }
}

/// How a job stands for one key and stamp.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Progress {
    Idle,
    Running,
    Done,
    Failed,
}

impl<K: Hash + Eq + Copy, S: PartialEq + Copy, T: Send + 'static> Jobs<K, S, T> {
    /// The result for `key` at `stamp`, or `None` while a worker computes
    /// it. When there is no result or job for this stamp, `prepare` runs
    /// on the calling thread to gather the worker's inputs and the worker
    /// it returns starts on a new thread; a job still running for an older
    /// stamp finishes first.
    pub fn get<W: FnOnce() -> T + Send + 'static>(
        &mut self,
        key: K,
        stamp: S,
        wake: Option<&egui::Context>,
        prepare: impl FnOnce() -> W,
    ) -> Option<Arc<T>> {
        self.poll(key);
        match self.slots.get(&key) {
            Some(Slot {
                stamp: s,
                state: State::Done(value),
            }) if *s == stamp => return Some(value.clone()),
            Some(Slot {
                stamp: s,
                state: State::Failed,
            }) if *s == stamp => return None,
            Some(Slot {
                state: State::Running(_),
                ..
            }) => return None,
            _ => {}
        }
        let (tx, rx) = channel();
        let wake = wake.cloned();
        let work = prepare();
        std::thread::spawn(move || {
            let _ = tx.send(work());
            if let Some(ctx) = wake {
                ctx.request_repaint();
            }
        });
        self.slots.insert(
            key,
            Slot {
                stamp,
                state: State::Running(rx),
            },
        );
        None
    }

    pub fn progress(&mut self, key: K, stamp: S) -> Progress {
        self.poll(key);
        match self.slots.get(&key) {
            Some(Slot { stamp: s, state }) if *s == stamp => match state {
                State::Running(_) => Progress::Running,
                State::Done(_) => Progress::Done,
                State::Failed => Progress::Failed,
            },
            Some(Slot {
                state: State::Running(_),
                ..
            }) => Progress::Running,
            _ => Progress::Idle,
        }
    }

    pub fn retain(&mut self, keep: impl Fn(&K) -> bool) {
        self.slots.retain(|k, _| keep(k));
    }

    fn poll(&mut self, key: K) {
        let Some(slot) = self.slots.get_mut(&key) else {
            return;
        };
        let State::Running(rx) = &slot.state else {
            return;
        };
        match rx.try_recv() {
            Ok(value) => slot.state = State::Done(Arc::new(value)),
            Err(TryRecvError::Disconnected) => slot.state = State::Failed,
            Err(TryRecvError::Empty) => {}
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::{Duration, Instant};

    fn wait(jobs: &mut Jobs<u32, u32, u32>, stamp: u32) -> Arc<u32> {
        let start = Instant::now();
        loop {
            if let Some(v) = jobs.get(1, stamp, None, || || 0) {
                return v;
            }
            assert!(
                start.elapsed() < Duration::from_secs(5),
                "job never finished"
            );
            std::thread::sleep(Duration::from_millis(2));
        }
    }

    #[test]
    fn a_result_arrives_later_and_is_kept_for_its_stamp() {
        let mut jobs: Jobs<u32, u32, u32> = Jobs::default();
        assert!(jobs.get(1, 7, None, || || 42).is_none());
        assert_eq!(*wait(&mut jobs, 7), 42);
        assert_eq!(jobs.progress(1, 7), Progress::Done);
        assert_eq!(jobs.progress(1, 8), Progress::Idle);
    }

    #[test]
    fn a_new_stamp_recomputes_after_the_old_job_finishes() {
        let mut jobs: Jobs<u32, u32, u32> = Jobs::default();
        assert!(jobs.get(1, 1, None, || || 10).is_none());
        let start = Instant::now();
        let value = loop {
            if let Some(v) = jobs.get(1, 2, None, || || 20) {
                break v;
            }
            assert!(start.elapsed() < Duration::from_secs(5));
            std::thread::sleep(Duration::from_millis(2));
        };
        assert_eq!(*value, 20);
    }

    #[test]
    fn a_panicking_worker_is_not_retried() {
        let mut jobs: Jobs<u32, u32, u32> = Jobs::default();
        assert!(jobs
            .get(1, 1, None, || || -> u32 { panic!("worker failed") })
            .is_none());
        let start = Instant::now();
        while jobs.progress(1, 1) == Progress::Running {
            assert!(start.elapsed() < Duration::from_secs(5));
            std::thread::sleep(Duration::from_millis(2));
        }
        assert_eq!(jobs.progress(1, 1), Progress::Failed);
        assert!(jobs.get(1, 1, None, || || 5).is_none());
        assert_eq!(jobs.progress(1, 1), Progress::Failed);
    }
}
