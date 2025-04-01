//! Sampling (oriented) point cloud surfaces.
use crate::sdf::SmoothDistanceField;
use nuts_rs::{CpuLogpFunc, LogpError};
use thiserror::Error;

#[derive(Debug, Error)]
#[error("Cannot query neighborhood for this sample")]
pub struct SurfaceSdfError;

impl LogpError for SurfaceSdfError {
    fn is_recoverable(&self) -> bool {
        true
    }
}

impl CpuLogpFunc for SmoothDistanceField {
    type LogpError = SurfaceSdfError;
    type TransformParams = ();

    fn dim(&self) -> usize {
        3
    }

    fn logp(&mut self, position: &[f64], grad: &mut [f64]) -> Result<f64, Self::LogpError> {
        let position: [f64; 3] = position.try_into().expect("position should be of length 3");

        let (smoothdist, gradient) = self
            .harmonic_gradient(position, self.surface_radius)
            .ok_or(SurfaceSdfError)?;

        grad.copy_from_slice(&gradient);

        // this is our log probability
        Ok(smoothdist)
    }
}
