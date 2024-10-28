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
    let hot_pink = Point3::new(1., 0., 0.5);
    let teal = Point3::new(0., 1., 0.5);
    let bspline = bspline_non_periodic();
    let (lower, upper) = (bspline.knots.lower, bspline.knots.upper);
    let low_v = bspline.knots.knots[lower].value;
    let upp_v = bspline.knots.knots[upper].value;
    let d = upp_v - low_v;

    let mut d0 = Vec::with_capacity(101);
    let mut d1 = Vec::with_capacity(101);
    let mut d2 = Vec::with_capacity(101);
    (0..=100)
        .for_each(|i| {
            let u = (i as f32) / 100f32 * d + low_v;
            let p = bspline.interop(u);
            let v1 = bspline.derivate_n(u, 1);
            let v2 = bspline.derivate_n(u, 2);
            let p0 = Point3::new(p.x, p.y, p.z);
            let p1 = Point3::new(v1.x, v1.y, v1.z + 0.5);
            let p2 = Point3::new(v2.x, v2.y, v2.z + 1.);

            d0.push(p0);
            d1.push(p1);
            d2.push(p2);
        });

    let axis = vec![
        (Point3::new(0., 0., 0.), Point3::new(1., 0., 0.)),
        (Point3::new(0., 0., 0.), Point3::new(0., 1., 0.)),
        (Point3::new(0., 0., 0.), Point3::new(0., 0., 1.)),
    ];

    while window.render() {
        axis.iter().for_each(|(p1, p2)| {
            window.draw_line(p1, p2, p2);
        });

        d0.iter().for_each(|pt| {
            window.draw_point(pt, &white);
        });
        d1.iter().for_each(|pt| {
            window.draw_point(pt, &hot_pink);
        });
        d2.iter().for_each(|pt| {
            window.draw_point(pt, &teal);
        });
    }
}

fn bspline_non_periodic() -> BSpline {
    let poles = vec![
        BSplinePole {
            pole: Vec3::new(0., 0., 0.,),
            weight: 2.
        },
        BSplinePole {
            pole: Vec3::new(0., 2., 0.,),
            weight: 1.
        },
        BSplinePole {
            pole: Vec3::new(1., 2., 0.,),
            weight: 3.
        },
        BSplinePole {
            pole: Vec3::new(3., 0., 0.,),
            weight: 1.
        },
        BSplinePole {
            pole: Vec3::new(4., 0., 0.,),
            weight: 2.
        },
        BSplinePole {
            pole: Vec3::new(4., 2., 0.,),
            weight: 4.
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
    BSpline {
        degree: 2,
        knots: BSplineKnots::try_new(knots, 2, poles.len(), false).unwrap(),
        poles,
        is_periodic: false,
        use_rational: true,
    }
}
