use super::bspline_error::{BSplineError, BSplineKnotsGenError};


#[derive(Debug, Clone)]
pub enum BSplineKnotsAlgo {
    DeBoor(Vec<f32>), // parameter
    Uniform(usize), // nb_poles
    Universal(usize), // nb_poles
}

#[derive(Debug, Clone, Copy)]
pub enum BSplineFeature {
    Periodic(usize), // front and last multiplicity times.
    Regular,
    ClampStart,
    ClampEnd,
    ClampAll,
}

#[derive(Debug, Clone, Copy)]
pub struct BSplineKnot {
    pub value: f32,
    pub multiplicity: usize,
}

#[derive(Debug, Clone)]
pub struct BSplineKnots {
    pub knots: Vec<BSplineKnot>,
    pub lower: usize,
    pub upper: usize,
    pub flatten: Vec<f32>,
    pub n_extend: usize,
}

impl BSplineKnots {
    pub fn check(knots: &[BSplineKnot], degree: usize, nb_poles: usize, is_periodic: bool)
    -> Result<(), BSplineError> {
        let m = knots.len();
        if m < 2 {
            return Err(BSplineError::TooFewKnots);
        }
        if nb_poles < 2 {
            return Err(BSplineError::TooFewPoles);
        }
        for i in 0..m {
            // check multiplicity
            let max_multiplicity = if (i == 0 || i == m - 1) && !is_periodic {
                degree + 1
            } else {
                degree
            };
            if !(knots[i].multiplicity <= max_multiplicity) {
                return Err(BSplineError::MultiplicityOverDegree);
            }
            // Knots should are increasing.
            if i > 0 {
                if knots[i].value <= knots[i - 1].value {
                    return Err(BSplineError::MultiplicityOverDegree);
                }
            }
        }
        if is_periodic {
            // First need to be the same with last.
            // First and last as same, then nb_knots should equal nb_poles.
            if knots[0].multiplicity != knots[m - 1].multiplicity {
                return Err(BSplineError::PeriodicEdgeNotMatch);
            }
            if nb_poles != knots.iter().skip(1).map(|k| k.multiplicity).sum::<usize>() {
                return Err(BSplineError::KnotsNumberNotMatch);
            }
        } else {
            // nb_knots = nb_poles + degree + 1.
            if nb_poles + degree + 1 != knots.iter().map(|k| k.multiplicity).sum::<usize>() {
                return Err(BSplineError::KnotsNumberNotMatch);
            }
        }
        Ok(())
    }

    fn non_periodic_bound(knots: &[BSplineKnot], degree: usize, at_start: bool) -> usize {
        let mut k_idx = if at_start { 0 } else { knots.len() - 1 };
        let mut n = 0;
        for _ in 0..=degree {
            n += 1;
            if n > knots[k_idx].multiplicity {
                n = 1;
                if at_start {
                    k_idx += 1;
                } else {
                    k_idx -= 1;
                }
                continue;
            }
        }
        k_idx
    }

    fn extend_knots(knots: &[BSplineKnot], degree: usize) -> Vec<f32> {
        let period = knots.last().unwrap().value - knots[0].value;
        let split = degree + 1 - knots[0].multiplicity;
        let mut extend = vec![0f32; 2 * split];
        let e_start = split - 1;
        let e_end = split;
        let mut i_start = knots.len() - 2;
        let mut i_end = 1;
        let mut n_start = 0;
        let mut n_end = 0;
        for i in 0..split {
            n_end += 1;
            if n_end > knots[i_end].multiplicity {
                i_end = if i_end == knots.len() - 1 {
                    0
                } else {
                    i_end + 1
                };
                n_end = 1;
            }
            extend[e_end + i] = knots[i_end].value + period;

            n_start += 1;
            if n_start > knots[i_start].multiplicity {
                i_start = if i_start == 0 {
                    knots.len() - 1
                } else {
                    i_start - 1
                };
                n_start = 1;
            }
            extend[e_start - i] = knots[i_start].value - period;
        }
        extend
    }

    pub fn try_new(
        knots: Vec<BSplineKnot>,
        degree: usize,
        nb_poles: usize,
        is_periodic: bool,
    ) -> Result<Self, BSplineError> {
        if let Err(e) = Self::check(&knots, degree, nb_poles, is_periodic) {
            return Err(e);
        }
        let (lower, upper) = if is_periodic {
            (0, knots.len() - 1)
        } else {
            let start = Self::non_periodic_bound(&knots, degree, true);
            let end = Self::non_periodic_bound(&knots, degree, false);
            (start, end)
        };
        let mut flatten = knots
            .iter()
            .flat_map(|knot| vec![knot.value; knot.multiplicity])
            .collect::<Vec<_>>();

        let n_extend = if is_periodic {
            let extend = Self::extend_knots(&knots, degree);
            let split = extend.len() / 2;
            let mut extended = vec![0f32; flatten.len() + extend.len()];
            extended[..split].clone_from_slice(&extend[..split]);
            extended[split..flatten.len() + split].clone_from_slice(&flatten);
            extended[flatten.len() + split..].clone_from_slice(&extend[split..]);
            flatten = extended;
            split
        } else { 0 };
        Ok(Self {
            knots,
            lower,
            upper,
            flatten,
            n_extend
        })
    }

    pub fn try_build(degree: usize, algo: &BSplineKnotsAlgo, feature: BSplineFeature) -> Result<BSplineKnots, BSplineKnotsGenError> {
        match algo {
            BSplineKnotsAlgo::DeBoor(vec) => new_de_boor(degree, vec, feature),
            BSplineKnotsAlgo::Uniform(nb_poles) => new_uniform(degree, *nb_poles, feature),
            BSplineKnotsAlgo::Universal(nb_poles) => new_universal(degree, *nb_poles, feature),
        }
    }
}

impl BSplineKnots {
    pub fn original_flatten(&self) -> Vec<f32> {
        self.knots.iter()
            .flat_map(|knot| {
                vec![knot.value; knot.multiplicity]
            }).collect()
    }

    pub fn get_knots_bounds(&self, knot_index: usize, degree: usize, is_periodic: bool) -> Vec<f32> {
        let idx = self
            .knots
            .iter()
            .take(knot_index + 1)
            .map(|k| k.multiplicity)
            .sum::<usize>();
        let idx = if is_periodic { idx } else {
            idx.max(degree).min(self.flatten.len()-degree-1)
        };
        let idx = idx + self.n_extend - 1;
        (idx - degree + 1..=idx + degree)
            .map(|i| self.flatten[i])
            .collect()
    }

    pub fn get_knot_index(&self, u: f32) -> usize {
        let mut idx = 0;
        for i in 1..self.knots.len() {
            if self.knots[i].value < u {
                idx += 1;
            }
        }
        idx.max(self.lower).min(self.upper)
    }

    pub fn get_flat_index(&self, u: f32) -> usize {
        let u_id = self.get_knot_index(u);
        let idx = self
            .knots
            .iter()
            .take(u_id + 1)
            .map(|k| k.multiplicity)
            .sum::<usize>() - 1;
        idx + self.n_extend
    }

    pub fn lower_value(&self) -> f32 {
        self.knots[self.lower].value
    }
    pub fn upper_value(&self) -> f32 {
        self.knots[self.upper].value
    }
}

fn new_uniform(degree: usize, nb_poles: usize, feature: BSplineFeature) -> Result<BSplineKnots, BSplineKnotsGenError> {
    let mut periodic = false;
    let mut total = degree + nb_poles + 1;
    let (first, last) = match feature {
        BSplineFeature::Periodic(s) => {
            periodic = true;
            total = 2 * nb_poles;
            (s, s)
        },
        BSplineFeature::Regular => (degree, degree),
        BSplineFeature::ClampStart => (degree+1, degree),
        BSplineFeature::ClampEnd => (degree, degree+1),
        BSplineFeature::ClampAll => (degree+1, degree+1),
    };
    if total < first + last {
        return Err(BSplineKnotsGenError::TooFewPoles);
    }
    let mid = total - first - last;
    let d = 1. / (mid + 1) as f32;

    let mut knots = vec![BSplineKnot { value: 0., multiplicity: first }];
    for i in 0..mid {
        knots.push(BSplineKnot { value: (i+1) as f32 * d, multiplicity: 1 });
    }
    knots.push(BSplineKnot { value: 1., multiplicity: last });

    match BSplineKnots::try_new(knots, degree, nb_poles, periodic) {
        Ok(knots) => Ok(knots),
        Err(e) => Err(BSplineKnotsGenError::BSplineErr(e)),
    }
}

fn new_universal(degree: usize, nb_poles: usize, feature: BSplineFeature) -> Result<BSplineKnots, BSplineKnotsGenError> {
    let mut periodic = false;
    let mut total = degree + nb_poles + 1;
    let (first, last) = match feature {
        BSplineFeature::Periodic(s) => {
            periodic = true;
            total = 2 * nb_poles;
            (s, s)
        },
        BSplineFeature::Regular => (degree, degree),
        BSplineFeature::ClampStart => (degree+1, degree),
        BSplineFeature::ClampEnd => (degree, degree+1),
        BSplineFeature::ClampAll => (degree+1, degree+1),
    };
    if total < first + last {
        return Err(BSplineKnotsGenError::TooFewPoles);
    }
    let mid = total - first - last;
    let factor = (nb_poles + 1 - degree) as f32;
    let mut knots = vec![BSplineKnot { value: 0., multiplicity: first }];
    for i in 0..mid {
        knots.push(BSplineKnot { value: (i+1) as f32 / factor, multiplicity: 1 });
    }
    knots.push(BSplineKnot { value: 1., multiplicity: last });

    match BSplineKnots::try_new(knots, degree, nb_poles, periodic) {
        Ok(knots) => Ok(knots),
        Err(e) => Err(BSplineKnotsGenError::BSplineErr(e)),
    }
}

pub fn new_de_boor(
    _degree: usize,
    _de_boor_param: &[f32],
    _feature: BSplineFeature
) -> Result<BSplineKnots, BSplineKnotsGenError>  {
    todo!()
}

// // nb_parameter = (n + 1)
// // nb_knots = (n + p + 1) = nb_parameter + degree
// pub fn new_de_boor_old(
//     degree: usize,
//     de_boor_param: Vec<f32>
// ) -> Option<Self> {
//     if de_boor_param.is_empty() {
//         return None;
//     }
//     if de_boor_param.len() == 1
//         && de_boor_param[0] >= 1. {
//         return None;
//     }
//     // check increasing
//     for i in 1..de_boor_param.len() {
//         if de_boor_param[i] < de_boor_param[i-1] {
//             return None;
//         }
//     }
//     let u0 = de_boor_param[0];
//     let u1 = *de_boor_param.last().unwrap();
//     let d = u1 - u0;
//     let de_boor_param = de_boor_param.into_iter().map(|v| v / d).collect::<Vec<f32>>();

//     let len = degree + de_boor_param.len();
//     let mut knots = vec![0f32; len];
//     let start = degree + 1;
//     let end = len - degree - 1;
//     for knot in end..len {
//         knots[knot] = d;
//     }
//     for k in start..end {
//         knots[k] = (k-degree..k-1).map(|i| de_boor_param[i]).sum::<f32>() / degree as f32 * d;
//     }
//     Some(Self { knots, start: u0, end: u1 })
// }