use super::bspline_poles::BSplinePole;
use f3l::glam::Vec3;

pub fn de_boor(
    u: f32,
    knots: &[f32],
    poles: &[BSplinePole],
    degree: usize,
    rational: bool,
) -> Vec3 {
    let mut rs = poles.iter().copied().collect::<Vec<_>>();

    for i in 0..degree {
        let mut rs_temp = vec![BSplinePole::default(); degree - i];
        for j in 0..degree - i {
            let a = (u - knots[j + i]) / (knots[j + degree] - knots[j + i]);
            rs_temp[j].pole = (1. - a) * rs[j].pole + a * rs[j + 1].pole;
            if rational {
                rs_temp[j].weight = (1. - a) * rs[j].weight + a * rs[j + 1].weight;
            }
        }
        rs = rs_temp;
    }
    if rational {
        rs[0].pole / rs[0].weight
    } else {
        rs[0].pole
    }
}

pub fn compute_coefficients(
    u: f32,
    knots: &[f32],
    nb_poles: usize,
    degree: usize,
    n_extend: usize,
) -> Vec<f32> {
    (0..nb_poles)
        .map(|i| recursive_coe(u, &knots, degree, i, n_extend))
        .collect()
}

fn recursive_coe(u: f32, knots: &[f32], p: usize, i: usize, n_extend: usize) -> f32 {
    if p == 0 {
        return if knots[n_extend + i] <= u && u < knots[n_extend + i + 1] {
            1f32
        } else {
            0f32
        };
    }

    let factor1 = knots[n_extend + i + p] - knots[n_extend + i];
    let factor2 = knots[n_extend + i + p + 1] - knots[n_extend + i + 1];
    let n_i_p_1 = recursive_coe(u, knots, p - 1, i, n_extend);
    let n_i_1_p_1 = recursive_coe(u, knots, p - 1, i + 1, n_extend);
    let part1 = if factor1 == 0. {
        0f32
    } else {
        (u - knots[n_extend + i]) * n_i_p_1 / factor1
    };
    let part2 = if factor2 == 0. {
        0f32
    } else {
        (knots[n_extend + i + p + 1] - u) * n_i_1_p_1 / factor2
    };
    part1 + part2
}
