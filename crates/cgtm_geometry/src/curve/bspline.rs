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

    pub fn derivate_n(&self, u:f32, d: usize) -> Vec3 {
        if d == 0 {
            return self.interop(u);
        }

        let mut knots = self.knots.clone();
        let mut poles = self.poles.clone();
        for deriv in 1..=d {
            let p = self.degree - deriv + 1;

            let flatten = knots.original_flatten();
            let mut temp_poles = vec![];
            for i in 0..poles.len()-1 {
                let factor = flatten[i+p+1] - flatten[i+1];
                let param = if factor==0f32 { 0f32 } else  { p as f32 / factor };
                temp_poles.push(BSplinePole {
                    pole: (poles[i+1].pole - poles[i].pole) * param,
                    weight: (poles[i+1].weight - poles[i].weight) * param
                });
            }
            poles = temp_poles;
            
            knots = knots.derivate(self.degree - deriv, self.poles.len()-deriv, self.is_periodic);
        }

        let bspline = BSpline {
            degree: self.degree - d,
            knots,
            poles,
            is_periodic: self.is_periodic,
            use_rational: self.use_rational,
        };
        bspline.interop(u)
    }
}

#[cfg(test)]
mod test_bspline {
    use bspline_knots::BSplineKnot;

    use super::*;

    #[test]
    fn derivate_1() {
        let poles = vec![
            BSplinePole {
                pole: Vec3::new(0., 0., 0.,),
                weight: 1.
            },
            BSplinePole {
                pole: Vec3::new(0., 2., 0.,),
                weight: 1.
            },
            BSplinePole {
                pole: Vec3::new(1., 2., 0.,),
                weight: 1.
            },
            BSplinePole {
                pole: Vec3::new(3., 0., 0.,),
                weight: 1.
            },
            BSplinePole {
                pole: Vec3::new(4., 0., 0.,),
                weight: 1.
            },
            BSplinePole {
                pole: Vec3::new(4., 2., 0.,),
                weight: 1.
            },
        ];
        let knots = vec![
            BSplineKnot {
                value: 0.,
                multiplicity: 3
            },
            BSplineKnot {
                value: 0.25,
                multiplicity: 1
            },
            BSplineKnot {
                value: 0.5,
                multiplicity: 1
            },
            BSplineKnot {
                value: 0.75,
                multiplicity: 1
            },
            BSplineKnot {
                value: 1.,
                multiplicity: 3
            },
        ];
        let bspline = BSpline {
            degree: 2,
            knots: BSplineKnots::try_new(knots, 2, poles.len(), false).unwrap(),
            poles,
            is_periodic: false,
            use_rational: false,
        };

        let lower = bspline.lower_parameter();
        let upper = bspline.upper_parameter();
        let d = upper - lower;
        let u = d * 0.5 + lower;

        let p = bspline.derivate_n(u, 1);
        assert_eq!(p, Vec3::new(8., -8., 0.));
        let p = bspline.derivate_n(u, 2);
        assert_eq!(p, Vec3::new(16., -32., 0.));
    }
}