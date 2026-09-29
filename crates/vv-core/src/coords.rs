//! One set of atom coordinates (a single model or trajectory frame).

use glam::Vec3;

#[derive(Clone, Debug, Default, PartialEq)]
pub struct CoordSet {
    positions: Vec<Vec3>,
}

impl CoordSet {
    pub fn new(positions: Vec<Vec3>) -> Self {
        Self { positions }
    }

    pub fn positions(&self) -> &[Vec3] {
        &self.positions
    }

    pub fn len(&self) -> usize {
        self.positions.len()
    }

    pub fn is_empty(&self) -> bool {
        self.positions.is_empty()
    }

    /// Axis-aligned bounding box as `(min, max)`, or `None` if empty.
    pub fn aabb(&self) -> Option<(Vec3, Vec3)> {
        let mut iter = self.positions.iter();
        let first = *iter.next()?;
        let (min, max) = iter.fold((first, first), |(lo, hi), p| (lo.min(*p), hi.max(*p)));
        Some((min, max))
    }

    /// Mean position: pulled toward denser regions, unlike
    /// the AABB midpoint `bounding_sphere` centres on -- the more
    /// representative "middle" of a molecule with an uneven atom count
    /// across its extent (a compact core with a floppy tail or branch).
    pub fn centroid(&self) -> Option<Vec3> {
        (!self.positions.is_empty())
            .then(|| self.positions.iter().sum::<Vec3>() / self.positions.len() as f32)
    }

    /// Bounding sphere centred on the AABB centre. Not minimal, but cheap and
    /// always enclosing — enough for camera framing and frustum tests.
    pub fn bounding_sphere(&self) -> Option<(Vec3, f32)> {
        let (min, max) = self.aabb()?;
        let center = (min + max) * 0.5;
        let radius_sq = self
            .positions
            .iter()
            .map(|p| p.distance_squared(center))
            .fold(0.0f32, f32::max);
        Some((center, radius_sq.sqrt()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_has_no_bounds() {
        let c = CoordSet::default();
        assert!(c.is_empty());
        assert_eq!(c.aabb(), None);
        assert_eq!(c.bounding_sphere(), None);
    }

    #[test]
    fn bounds_enclose_all_points() {
        let c = CoordSet::new(vec![
            Vec3::new(-1.0, 0.0, 0.0),
            Vec3::new(3.0, 2.0, -5.0),
            Vec3::new(0.0, 4.0, 1.0),
        ]);
        let (min, max) = c.aabb().unwrap();
        assert_eq!(min, Vec3::new(-1.0, 0.0, -5.0));
        assert_eq!(max, Vec3::new(3.0, 4.0, 1.0));
        let (center, radius) = c.bounding_sphere().unwrap();
        assert_eq!(center, Vec3::new(1.0, 2.0, -2.0));
        for p in c.positions() {
            assert!(p.distance(center) <= radius + 1e-6);
        }
    }
}
