use bspline_knots::BSplineKnots;
use bspline_poles::BSplinePole;
use f3l::glam::Vec3;
use utils::{compute_coefficients, de_boor};

pub mod bspline_approximation;
pub mod bspline_data_knots;
pub mod bspline_error;
pub mod bspline_knots;
pub mod bspline_poles;
pub mod utils;

#[derive(Debug, Clone)]
pub struct BSpline {
    pub degree: usize,
    pub knots: BSplineKnots,
    pub poles: Vec<BSplinePole>,
    pub is_periodic: bool,
    pub use_rational: bool,
}

impl BSpline {
    fn find_pole_index(
        knots: &BSplineKnots,
        knot_index: usize,
        degree: usize,
        is_periodic: bool,
    ) -> usize {
        let pole_index = knots
            .knots
            .iter()
            .take(knot_index + 1)
            .map(|knot| knot.multiplicity)
            .sum::<usize>();
        if is_periodic {
            pole_index - knots.knots[knots.lower].multiplicity
        } else {
            pole_index - (degree + 1)
        }
    }

    fn get_poles(&self, pole_index: usize) -> Vec<BSplinePole> {
        let mut pi = pole_index;
        (0..=self.degree)
            .map(|_| {
                if pi > self.poles.len() - 1 {
                    pi = 0;
                }
                let mut p = self.poles[pi];
                if self.use_rational {
                    p.pole *= p.weight;
                }
                pi += 1;
                p
            })
            .collect()
    }

    pub fn interop(&self, u: f32) -> Vec3 {
        let knot_index = self.knots.get_knot_index(u);
        let knot_slice = self
            .knots
            .get_knots_bounds(knot_index, self.degree, self.is_periodic);
        let pole_index =
            Self::find_pole_index(&self.knots, knot_index, self.degree, self.is_periodic);
        let poles = self.get_poles(pole_index);

        de_boor(u, &knot_slice, &poles, self.degree, self.use_rational)
    }

    pub fn lower_parameter(&self) -> f32 {
        self.knots.lower_value()
    }

    pub fn upper_parameter(&self) -> f32 {
        self.knots.upper_value()
    }

    pub fn coefficients(&self, u: f32) -> Vec<f32> {
        compute_coefficients(
            u,
            &self.knots.flatten,
            self.poles.len(),
            self.degree,
            self.knots.n_extend,
        )
    }
}
