use f3l::glam::Vec3;
use na::{Dyn, OMatrix, RowVector3, U3};

use crate::curve::bspline::bspline_knots::BSplineKnots;

use super::{bspline_data_knots::BSplineDataKnotsAlgo, bspline_error::BSplineApproximationError, bspline_knots::{BSplineFeature, BSplineKnotsAlgo}, bspline_poles::BSplinePole, utils::compute_coefficients, BSpline};


#[derive(Debug, Clone)]
pub struct BSplineApproximation {
    pub degree: usize,
    pub nb_poles: usize,
    pub feature: BSplineFeature,
    pub knots_algo: BSplineKnotsAlgo,
    pub data_knots_algo: BSplineDataKnotsAlgo,
}

impl BSplineApproximation {
    pub fn try_approximate(&self, data: &[Vec3]) -> Result<BSpline, BSplineApproximationError> {
        if let BSplineFeature::Periodic(_) = self.feature {
            return Err(BSplineApproximationError::NotAllowPeriodic);
        }

        let bspline_knots = match BSplineKnots::try_build(
            self.degree,
            &self.knots_algo,
            self.feature
        ) {
            Ok(knots) => knots,
            Err(e) => return  Err(BSplineApproximationError::BSplineKnotsGeneration(e)),
        };
        let knots = bspline_knots.original_flatten();
        let u_vector = self.data_knots_algo.generate(data);

        let h = self.nb_poles - 1;
        let n = data.len() - 1;

        let ns = u_vector
            .iter()
            .map(|&u| {
                // only allow non-periodic feature, so n_extend must be 0.
                compute_coefficients(u, &knots, self.nb_poles, self.degree, 0)
            })
            .collect::<Vec<_>>();

        let mut qks = (1..=n - 1)
            .map(|k| {
                let coe = &ns[k];
                data[k] - coe[0] * data[0] - coe[h] * data[n]
            })
            .collect::<Vec<_>>();
        qks.insert(0, Vec3::ZERO);

        let q_rows = (1..=h - 1)
            .map(|i| {
                let r = (1..=n - 1).map(|k| ns[k][i] * qks[k]).sum::<Vec3>();
                RowVector3::new(r.x, r.y, r.z)
            })
            .collect::<Vec<_>>();
        let q_matrix = OMatrix::<f32, Dyn, U3>::from_rows(&q_rows);
        let mut n_matrix = OMatrix::<f32, Dyn, Dyn>::zeros(n - 1, h - 1);
        (1..=n - 1).for_each(|k| {
            (1..=h - 1).for_each(|i| *n_matrix.index_mut((k - 1, i - 1)) = ns[k][i]);
        });
        let m_matrix = n_matrix.transpose() * n_matrix;
        let decompose = m_matrix.svd(true, true);
        let p_matrix = decompose.solve(&q_matrix, 1e-6).expect("Lu Solve (Ax=B).");

        let mut poles = vec![data[0]];
        for row in p_matrix.row_iter() {
            poles.push(Vec3::new(row[0], row[1], row[2]));
        }
        poles.push(data[n]);

        let bspline_poles = poles.into_iter()
            .map(|pole| {
                BSplinePole {pole, weight: 1f32}
            }).collect::<Vec<_>>();
        
        Ok(BSpline {
            degree: self.degree,
            is_periodic: false,
            knots: bspline_knots,
            poles: bspline_poles,
            use_rational: false
        })
    }
}