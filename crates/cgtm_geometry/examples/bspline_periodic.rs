use cgtm_geometry::curve::bspline::bspline_knots::{BSplineKnot, BSplineKnots};
use cgtm_geometry::curve::bspline::bspline_poles::BSplinePole;
use cgtm_geometry::curve::bspline::BSpline;
use f3l::glam::Vec3;
#[cfg(feature = "app")]
use kiss3d::light::Light;
#[cfg(feature = "app")]
use kiss3d::window::Window;

#[cfg(feature = "app")]
fn main() {
    use kiss3d::nalgebra::Point3;

    let mut window = Window::new("Kiss3d: points");

    window.set_light(Light::StickToCamera);
    window.set_point_size(10.0);

    let white = Point3::new(1., 1., 1.);
    let bspline = bspline_periodic();
    let (lower, upper) = (bspline.knots.lower, bspline.knots.upper);
    let low_v = bspline.knots.knots[lower].value;
    let upp_v = bspline.knots.knots[upper].value;
    let d = upp_v - low_v;
    let pts = (0..=100)
        .map(|i| {
            let u = (i as f32) / 100f32 * d + low_v;
            let p = bspline.interop(u);
            Point3::new(p.x, p.y, p.z)
        })
        .collect::<Vec<_>>();

    while window.render() {
        pts.iter().for_each(|pt| {
            window.draw_point(pt, &white);
        });
    }
}

fn bspline_periodic() -> BSpline {
    let degree = 3_usize;
    let nb_poles = 9_usize;
    let is_periodic = true;
    let knots = vec![
        BSplineKnot {
            value: 0.,
            multiplicity: 3,
        },
        BSplineKnot {
            value: 1.,
            multiplicity: 2,
        },
        BSplineKnot {
            value: 3.,
            multiplicity: 2,
        },
        BSplineKnot {
            value: 6.,
            multiplicity: 2,
        },
        BSplineKnot {
            value: 10.,
            multiplicity: 3,
        },
    ];
    let knots = BSplineKnots::try_new(knots, degree, nb_poles, is_periodic);
    let knots = knots.unwrap();

    let poles = vec![
        BSplinePole {
            pole: Vec3::new(1., 0., 0.),
            weight: 1.,
        },
        BSplinePole {
            pole: Vec3::new(1., 1., 0.),
            weight: 1.,
        },
        BSplinePole {
            pole: Vec3::new(0., 1., 0.),
            weight: 1.,
        },
        BSplinePole {
            pole: Vec3::new(-1., 1., 0.),
            weight: 1.,
        },
        BSplinePole {
            pole: Vec3::new(-1., 0., 0.),
            weight: 1.,
        },
        BSplinePole {
            pole: Vec3::new(-1., -1., 0.),
            weight: 1.,
        },
        BSplinePole {
            pole: Vec3::new(0., -1., 0.),
            weight: 1.,
        },
        BSplinePole {
            pole: Vec3::new(1., -1., 0.),
            weight: 1.,
        },
        BSplinePole {
            pole: Vec3::new(1., 0., 0.),
            weight: 1.,
        },
    ];

    BSpline {
        degree,
        is_periodic: true,
        knots,
        poles,
        use_rational: false,
    }
}
