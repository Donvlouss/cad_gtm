use f3l::glam::Vec3;

#[derive(Debug, Clone, Copy, Default)]
pub enum BSplineDataKnotsAlgo {
    #[default]
    Uniform,
    ChordLength,
    Centripetal(f32),
}

impl BSplineDataKnotsAlgo {
    pub const CENTRIPETAL_SQRT: f32 = 0.5;

    pub fn generate(&self, data: &[Vec3]) -> Vec<f32> {
        match self {
            BSplineDataKnotsAlgo::Uniform => uniform(data),
            BSplineDataKnotsAlgo::ChordLength => centripetal(data, 1f32),
            BSplineDataKnotsAlgo::Centripetal(alpha) => centripetal(data, *alpha),
        }
    }
}

fn uniform(data: &[Vec3]) -> Vec<f32> {
    let n = 1f32 / data.len() as f32;
    (0..data.len()).map(|i| (i + 1) as f32 * n).collect()
}

fn centripetal(data: &[Vec3], alpha: f32) -> Vec<f32> {
    let (mut acc, sum) = data.iter().enumerate().fold(
        (Vec::<f32>::with_capacity(data.len()), 0f32),
        |(mut acc, mut sum), (i, &v)| {
            if i == 0 {
                acc.push(0f32);
            } else {
                let d = (v.distance(data[i - 1])).powf(alpha);
                acc.push(acc.last().unwrap() + d);
                sum += d;
            }
            (acc, sum)
        },
    );
    for a in acc.iter_mut() {
        *a /= sum;
    }
    acc
}
