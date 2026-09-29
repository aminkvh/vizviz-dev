//! Adaptive quality: trades impostor detail for frame time under load.

/// Controls the projected-radius threshold (pixels) below which an atom is
/// drawn as a point instead of a ray-cast sphere. Raising it moves atoms
/// from the expensive bucket to the cheap one. The cull shader jitters the
/// threshold per atom, so the quad/point split changes smoothly as this
/// value moves instead of flipping all at once.
#[derive(Clone, Debug)]
pub struct AdaptiveLod {
    pub threshold_px: f32,
    pub min_px: f32,
    pub max_px: f32,
    /// Frame time the controller steers towards.
    pub target_ms: f32,
    /// Below this the threshold is allowed to relax (more detail).
    pub relax_ms: f32,
}

impl Default for AdaptiveLod {
    fn default() -> Self {
        Self {
            threshold_px: 1.0,
            min_px: 1.0,
            max_px: 24.0,
            target_ms: 25.0,
            relax_ms: 15.0,
        }
    }
}

impl AdaptiveLod {
    /// Feed the last frame's GPU time; returns the threshold to use next.
    /// Proportional: a frame twice the target raises the threshold by ~20%,
    /// a frame at half the relax point lowers it by ~10%.
    pub fn update(&mut self, frame_ms: f32) -> f32 {
        let ratio = if frame_ms > self.target_ms {
            (frame_ms / self.target_ms).sqrt().min(1.2)
        } else if frame_ms < self.relax_ms {
            (frame_ms / self.relax_ms).sqrt().max(0.9)
        } else {
            1.0
        };
        self.threshold_px = (self.threshold_px * ratio).clamp(self.min_px, self.max_px);
        self.threshold_px
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn raises_under_load_and_relaxes_when_fast() {
        let mut lod = AdaptiveLod::default();
        let slow = lod.update(40.0);
        assert!(slow > 1.0 && slow <= 1.2);
        for _ in 0..100 {
            lod.update(100.0);
        }
        assert_eq!(lod.threshold_px, lod.max_px);
        for _ in 0..200 {
            lod.update(5.0);
        }
        assert_eq!(lod.threshold_px, lod.min_px);
    }

    #[test]
    fn holds_inside_the_deadband() {
        let mut lod = AdaptiveLod {
            threshold_px: 4.0,
            ..Default::default()
        };
        assert_eq!(lod.update(20.0), 4.0);
        assert_eq!(lod.update(15.0), 4.0);
        assert_eq!(lod.update(25.0), 4.0);
    }

    #[test]
    fn steps_are_proportional_and_bounded() {
        let mut lod = AdaptiveLod {
            threshold_px: 4.0,
            ..Default::default()
        };
        let mild = lod.update(30.0);
        assert!(mild > 4.0 && mild < 4.5);
        let mut lod = AdaptiveLod {
            threshold_px: 4.0,
            ..Default::default()
        };
        let severe = lod.update(400.0);
        assert_eq!(severe, 4.8);
    }
}
