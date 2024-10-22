use cgtm_geometry::curve::bspline::bspline_approximation::BSplineApproximation;
use cgtm_geometry::curve::bspline::bspline_data_knots::BSplineDataKnotsAlgo;
use cgtm_geometry::curve::bspline::bspline_knots::{BSplineFeature, BSplineKnotsAlgo};
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

    let white = Point3::new(1., 1., 1.,);
    let green = Point3::new(0., 1., 0.);
    let data = gen_data();
    let bspline = approximation(&data);
    let (lower, upper) = (bspline.knots.lower, bspline.knots.upper);
    let low_v = bspline.knots.knots[lower].value;
    let upp_v = bspline.knots.knots[upper].value;
    let d = upp_v - low_v;

    let ori = data.iter().map(|p| Point3::new(p.x, p.y, p.z)).collect::<Vec<_>>();
    let pts = (0..400)
        .map(|i| {
            let u = (i as f32) / 400f32 * d + low_v;
            let p = bspline.interop(u);
            Point3::new(p.x, p.y, p.z)
        }).collect::<Vec<_>>();

    while window.render() {
        pts.iter().for_each(|pt| {
            window.draw_point(pt, &green);
        });
        ori.iter().for_each(|pt| {
            window.draw_point(pt, &white);
        });
    }
}

fn gen_data() -> Vec<Vec3> {
    let samples = 400;
    let start = -2f32;
    let end = 2f32;
    let d = (end - start) / samples as f32;
    (0..samples).map(|x| {
        let x = x as f32 * d;
        let y = x.sin();
        Vec3::new(x, y, 0.)
    }).collect::<Vec<_>>()
}

fn approximation(data: &[Vec3]) -> BSpline {

    let approx = BSplineApproximation {
        degree: 2,
        nb_poles: 9,
        feature: BSplineFeature::ClampAll,
        knots_algo: BSplineKnotsAlgo::Uniform(9),
        data_knots_algo: BSplineDataKnotsAlgo::Uniform,
    };
    match approx.try_approximate(data) {
        Ok(b) => b,
        Err(e) => panic!("{e}"),
    }
}