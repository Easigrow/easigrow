//! Calculate various plastic zone properties.

use std::f64::consts::PI;
use std::ops::{Add, Sub};

#[derive(Debug, Clone)]
pub struct ZoneWidth {
    pub plane_stress: f64,
    pub plane_strain: f64,
}

impl Add<f64> for ZoneWidth {
    type Output = ZoneWidth;

    fn add(self, offset: f64) -> ZoneWidth {
        ZoneWidth {
            plane_stress: self.plane_stress + offset,
            plane_strain: self.plane_strain + offset,
        }
    }
}

impl Sub<f64> for ZoneWidth {
    type Output = ZoneWidth;

    fn sub(self, offset: f64) -> ZoneWidth {
        ZoneWidth {
            plane_stress: self.plane_stress - offset,
            plane_strain: self.plane_strain - offset,
        }
    }
}

/// Calculate the size of the plastic zone for (plane stress, plane strain)
/// Ref. Anderson P. 485
/// Ref. Tada 1973 P. 1.17
pub fn zone_size(kmax: f64, sigma_yield: f64) -> ZoneWidth {
    // alpha is a constraint factor
    let alpha_plane_stress = 2.0;
    let alpha_plane_strain = 6.0;

    // radius of onset of yield with adjustment of the shape of the yield zone
    let radius_yield = |alpha: f64| (1.0 / (alpha * PI)) * (kmax / sigma_yield).powi(2);

    // the plastic zone width is approximately 2 * radius_yield
    ZoneWidth {
        plane_stress: 2.0 * radius_yield(alpha_plane_stress),
        plane_strain: 2.0 * radius_yield(alpha_plane_strain),
    }
}

pub fn zone_size_new(kmax: f64, sigma_yield: f64, alpha: f64) -> f64 {
    (1.0 / PI) * (kmax / (alpha * sigma_yield)).powi(2)
}

/// Use an iterative approach to calculate plastic zone length in configurations
/// with a hole. Configurations without a hole will skip the iterations.
// Note that the iterative code for a hole has not been implemented
pub fn zone_size_iterative(
    s_max: f64,
    crack_length: f64,
    spectrum_max: f64,
    plate_width: f64,
    flow_stress: f64,
) -> Result<f64, String> {
    // Don't bother if s_max is relatively small
    if s_max < 0.05 * spectrum_max {
        return Ok(0.0);
    }

    if crack_length > plate_width {
        return Err("Crack length was greater than plate width.".to_string());
    }

    let picw = 0.5 * PI * crack_length / plate_width;
    let snc = picw.sin();

    let scs = snc / (PI * s_max / (2.0 * flow_stress)).cos();

    if scs >= 1.0 {
        return Err("Plastic zone bounds error.".to_string());
    }

    let r = crack_length * (scs.asin() / picw - 1.0);
    let d = crack_length + r;

    if d >= plate_width {
        return Err("Crack length + plastic zone was greater than plate width.".to_string());
    }

    Ok(r)
}

#[cfg(test)]
mod tests {
    use super::zone_size_iterative;

    #[test]
    fn zone_size_iterative_returns_zero_when_s_max_is_small() {
        let s_max = 4.9;
        let spectrum_max = 100.0;

        // The below values should have no effect and are set arbitrarily
        let crack_length = 1.0e-4;
        let plate_width = 1.0e-2;
        let flow_stress = 400.0;

        let result =
            zone_size_iterative(s_max, crack_length, spectrum_max, plate_width, flow_stress);

        assert_eq!(result, Ok(0.0));
    }

    #[test]
    fn zone_size_iterative_returns_error_when_crack_length_greater_than_plate_width() {
        let crack_length = 2.0e-2;
        let plate_width = 1.0e-2;

        // s_max needs to be >= 5% of spectrum max
        let s_max = 90.0;
        let spectrum_max = 100.0;

        // The below values should have no effect and are set arbitrarily
        let flow_stress = 400.0;

        let result =
            zone_size_iterative(s_max, crack_length, spectrum_max, plate_width, flow_stress);

        assert!(result.is_err());
    }

    #[test]
    fn zone_size_iterative_is_correct_for_normal_inputs() {
        let crack_length = 1.0e-4;
        let plate_width = 0.0125;
        let s_max = 250.0;
        let spectrum_max = 300.0;
        let flow_stress = 400.0;

        let result =
            zone_size_iterative(s_max, crack_length, spectrum_max, plate_width, flow_stress).unwrap();
        assert!((result - 8.00058e-5).abs() < 1e-10);
    }
}
