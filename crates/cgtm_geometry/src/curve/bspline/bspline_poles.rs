use f3l::glam::Vec3;

#[derive(Debug, Clone, Copy, Default)]
pub struct BSplinePole {
    pub pole: Vec3,
    pub weight: f32,
}