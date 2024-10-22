use thiserror::Error;

#[derive(Debug, Clone, Copy, Error)]
pub enum BSplineError {
    #[error("Knots should be greater than 2.")]
    TooFewKnots,
    #[error("Poles should be greater than 2.")]
    TooFewPoles,
    #[error("Multiplicity should be less or equal `degree`. If non-periodic, first and last could be `degree+1`.")]
    MultiplicityOverDegree,
    #[error("When periodic, multiplicity of first and last should be equal.")]
    PeriodicEdgeNotMatch,
    #[error("When non-periodic, knots number should be poles number + degree + 1. When periodic, knots number - last multiplicity = poles number.")]
    KnotsNumberNotMatch,
}

#[derive(Debug, Clone, Copy, Error)]
pub enum BSplineKnotsGenError {
    #[error(transparent)]
    BSplineErr(BSplineError),
    #[error("Too few poles to generate knots.")]
    TooFewPoles,
}

#[derive(Debug, Clone, Copy, Error)]
pub enum BSplineApproximationError {
    #[error(transparent)]
    BSplineKnotsGeneration(#[from] BSplineKnotsGenError),
    #[error("Not allow to approximate periodic bspline.")]
    NotAllowPeriodic,
}
