//! This module optimises the model to best match crack growth measurements.
//!

extern crate rayon;

use self::rayon::prelude::*;
use fatigue::grow::{FatigueTest, BasicHistory};
use log::debug;
use std::f64;

use fatigue::io::Measurement;
use fatigue::grow;

pub mod nelder;
pub mod particle_swarm;
pub mod sweep;

pub trait Optimisation {
    fn run(&mut self, samples: &mut [Sample], use_log_crack_depth: bool);

    fn get_parameters(&self) -> &[f64];

    fn get_result(&self) -> String;
}

/// A single optimisation sample. Some optimisation runs
/// may include many samples.
#[derive(Debug, Clone)]
pub struct Sample {
    /// Predicted growth
    pub fatigue_test: Box<dyn FatigueTest + Send + Sync>,
    /// Fractography measurements to be matched
    pub measurements: Vec<Measurement>,
    /// Change the impact this test has on the overall result
    pub weight: f64,
}

/// Calculates the total error between predictions and measurements.
pub fn sum_prediction_error_factorised(samples: &mut [Sample], params: &[f64], factors: &[f64], use_log_crack_depth: bool) -> f64 {
    let mut floored_params = Vec::new();
    // factors is a nondimensionalised parameter
    for (i, parameter) in params.iter().enumerate() {
        floored_params.push(parameter * factors[i].max(0.0)); // restrict variable from going negative - (why? CH 22-11-11)
    }

    sum_prediction_error(samples, &floored_params, use_log_crack_depth)
}

/// Calculates the total error between predictions and measurements.
pub fn _sum_prediction_error_normalised(samples: &mut [Sample], params: &[f64], scale: f64, use_log_crack_depth: bool) -> f64 {
    let mut floored_params = Vec::new();
    // factors is a nondimensionalised parameter
    for (i, _parameter) in params.iter().enumerate() {
        floored_params.push(params[i].max(0.0) * scale); // restrict variable from going negative - (why? CH 22-11-11)
    }

    sum_prediction_error(samples, &floored_params, use_log_crack_depth)
}

/// Calculates the total error between predictions and measurements.
pub fn sum_prediction_error(samples: &mut [Sample], params: &[f64], use_log_crack_depth: bool) -> f64 {
    let all_errors = prediction_error(params, samples, use_log_crack_depth);
    let mut result = 0.0;

    for (errors, sample) in all_errors.iter().zip(samples) {
        let mut sum = 0.0;
        for error in errors {
            sum += error * error;
        }
        result += sum * sample.weight
    }
    result.sqrt()
}

/// Error function is the difference between predicted and measured crack growth rates.
///
/// This is typically the objective function to be minimised.
pub fn prediction_error(params: &[f64], samples: &mut [Sample], use_log_crack_depth: bool) -> Vec<Vec<f64>> {
    let mut errors: Vec<Vec<f64>> = Vec::new();

    // calculate the errors
    debug!("Using the material parameters {:?}", params);

    samples
        .par_iter_mut()
        .map(|sample| {
            sample.fatigue_test.update_parameters(params);
            let history = sample.fatigue_test.run_for_optimisation();

            let result = match compare_growth(&sample.measurements, history, use_log_crack_depth) {
                Some(result) => result,
                None => vec![f64::INFINITY, 1.0],
            };

            // We need to clear the cached history, otherwise we'll run into memory issues
            sample.fatigue_test.reset();

            result
        })
        .collect_into_vec(&mut errors);

    let mut total_error = Vec::new();
    for (_sample, prediction_error) in samples.iter().zip(&errors) {
        let rms_error = rms(prediction_error);
        // println!("sample error: {}", rms_error);
        total_error.push(rms_error);
    }

    debug!("Total error: {}", rms(&total_error));

    errors
}

/// Calculate the Root-Mean-Square value of a sequence.
pub fn rms(x: &[f64]) -> f64 {
    x.iter().fold(0.0f64, |acc, &a| acc + a.powi(2)).sqrt()
}

/// Normalise a list of numbers. Returns the scale.
pub fn _normalise(x: &mut [f64]) -> f64 {
    let scale = (*x.iter().max_by(|a, b| a.abs().partial_cmp(&b.abs()).unwrap()).unwrap()).abs();

    for val in x {
        *val /= scale;
    }

    scale
}

/// Calculate the errors between predicted and measured data.
fn compare_growth(
    measured: &[Measurement],
    history: &[grow::BasicHistory],
    use_log_crack_depth: bool,
) -> Option<Vec<f64>> {
    // Not enough data to compare against, most likely reason is the dadn parameters are
    // very inaccurate
    if history.len() <= 1 {
        return None;
    }

    let mut errors = Vec::with_capacity(measured.len());

    // Find the predicted crack growth rates at the measured crack
    // size by looking for the corresponding history that match the
    // two adjacent fracto measurements.
    for i in 0..(measured.len() - 1) {
        // Skip the comparison if crack length is zero implying a
        // discontinuous measurement.
        if (measured[i].a == 0.0) || (measured[i + 1].a == 0.0) {
            continue;
        }

        let h_start = find_nearest(measured[i].a, history);
        let h_end = find_nearest(measured[i + 1].a, history);

        let history_dadb = if h_start == h_end {
            0.0
        } else {
            let mut h_end_value;
            let mut h_start_value;
            h_end_value = history[h_end].crack_length;
            h_start_value = history[h_start].crack_length;
            
            if use_log_crack_depth {
                h_end_value = h_end_value.ln();
                h_start_value = h_start_value.ln();
            }
            
            // da / dB
            (h_end_value - h_start_value) / (history[h_end].block - history[h_start].block)
        };

        let measured_da = if use_log_crack_depth {
            measured[i + 1].a.ln() - measured[i].a.ln()
        } else {
            measured[i + 1].a - measured[i].a
        };
        
        let measured_dadb = measured_da / (measured[i + 1].block - measured[i].block);

        // if (history[h_start].geometry.a.as_ref().unwrap().length - 0.000457882221059177).abs() <= f64::EPSILON {
        // if (history[h_end].block - 21.0).abs() < 0.0002 {
        //     for j in h_end - 3..=h_end + 2 {
        //         println!("history[{}]: {}", j, history[j].geometry.a.as_ref().unwrap().length);
        //     }

        //     println!("measured[{}].a:  {}", i + 1, measured[i + 1].a);

        //     println!("history[h_end].block: {}, history[h_start].block: {}", history[h_end].block, history[h_start].block);
        //     println!("history[h_end].block: {}, history[h_start].block: {}", (history[h_end].block * 10000.0).round() / 10000.0, (history[h_start].block * 10000.0).round() / 10000.0);
        //     println!("measured[i + 1].block: {}, measured[i].block: {}", measured[i + 1].block, measured[i].block);
        //     println!("history_da: {}", history[h_end].crack_length - history[h_start].crack_length);
        //     println!("history[h_end - 1].a.length: {}, block: {}", history[h_end - 1].crack_length, history[h_end - 1].block );
        //     println!("history[h_end].a.length: {}", history[h_end].crack_length );
        //     println!("history[h_start].a.length: {}", history[h_start].crack_length );
        //     println!("history[h_start + 1].a.length: {}", history[h_start + 1].crack_length );
        //     println!("measured_da: {}", measured_da);
        //     println!("history_dadb: {}", history_dadb);
        //     println!("measured_dadb: {}", measured_dadb);
        //     println!("hend: {}, hstart: {}, his.len(): {}", h_end, h_start, history.len());
        // }
        

        // The choice of error term is difficult. Here we use the
        // difference between logs, which says that the order of
        // magnitude of the accuracy is important with equal
        // weighting being given to being out by an order of
        // magnitude at high crack growth rates as with small
        // crack growth rates. In the absence of any knowledge on
        // the number of cycles of each this seems like the best
        // error function.
        let mut error = history_dadb.ln() - measured_dadb.ln();

        // It's possible for this to occur in some cases where a dadn returns a
        // negative value as a result of being provided bad parameters from the
        // optimiser. This can cause a negative history delta, which then fails
        // if we try to take the log.
        if error.is_nan() || error.is_infinite() || history_dadb == 0.0 {
            error = -100.0;
        }

        // In addition to growth rate, we will compare absolute distance
        // between the two curves. It's possible for a prediction to diverge
        // early and then after a few blocks the growth rates may become similar.
        // This results in curves that, to the eye, are off by a lot. By taking
        // absolute distance between the measured and predicted, we can penalise
        // these results.

        // Need to subtract the initial block value from the measurements, as it may not be 0
        // let start_block_error = (history[h_start].block - measured[i].block - measured[0].block).abs();
        // let end_block_error = (history[h_end].block - measured[i + 1].block - measured[0].block).abs();

        // Early deviations have a higher impact on the result
        // let position_bias = 1.0 - i as f64 / (measured.len() - 1) as f64;

        // error += (start_block_error + end_block_error) * position_bias;

        // if *verbosity == options::Verbosity::Verbose {
        //     println!("New Measurement comparison summary:");
        //     println!(
        //         "  History: start {:6.2} at a={}, end {:6.2} at a={}",
        //         history[h_start].block,
        //         history[h_start].geometry.c.as_ref().unwrap().length,
        //         history[h_end].block,
        //         history[h_end].geometry.c.as_ref().unwrap().length
        //     );
        //     println!(
        //         "  Measurement: start {} at a={}, end {} at a={}",
        //         measured[i].block,
        //         measured[i].a,
        //         measured[i + 1].block,
        //         measured[i + 1].a
        //     );
        //     println!(
        //         "  giving average growth rates: history {:8e}, measured {:8e}, error {:8e}",
        //         history_dadb, measured_dadb, error
        //     );
        // }

        errors.push(error);
    }

    Some(errors)
}

/// Find the block in history which has a crack length closest to the target
/// Returns the index of the closest match
fn find_nearest(target: f64, history: &[BasicHistory]) -> usize {

    // Handle the case that the crack grows "infinitely" within the first block
    if history[0].crack_length.is_infinite() {
        return 0;
    }
    
    let comparison_function = 
        |probe: &BasicHistory| probe.crack_length.partial_cmp(&target).unwrap();
    
    let mut result = match history.binary_search_by(comparison_function) {
        // Found an exact match, return the exact value
        Ok(index) => {
            index
        },
        Err(index) => {
            if index == history.len() {
                index - 1
            } else if index > 0 {
                let lower = (target - history[index - 1].crack_length).abs();
                let upper = (target - history[index].crack_length).abs();

                if lower < upper {
                    index - 1
                } else {
                    index
                }
            } else {
                index
            }
        }
    };

    // Binary search will return the index where an element could be inserted
    // and maintain order. In some cases, this index will be after the last
    // item. This function must only return valid indexes, so in this case
    // return the index of the last item.
    if result == history.len() {
        result -= 1;
    }

    result
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeMap;
    use std::f64::EPSILON;
    use std::f64::consts::FRAC_PI_2;

    use crate::options;

    use super::{compare_growth, _normalise, sum_prediction_error};
    use fatigue::beta::BetaResult;
    use fatigue::dadn::ParameterLabel;
    use fatigue::grow::{CrackLengthDelta, CrackGeometry, CrackFront, StressIntensity, StressIntensityDelta};
    use fatigue::io::Measurement;
    use fatigue::{cycle, grow, tag};

    #[test]
    fn test_compare_growth() {
        let measured = vec![
            Measurement {
                line: 0,
                block: 20.0,
                a: 0.11,
            },
            Measurement {
                line: 0,
                block: 30.0,
                a: 0.20,
            },
            Measurement {
                line: 0,
                block: 40.0,
                a: 0.31,
            },
        ];

        let history = vec![
            fake_history(110.0, 0.05),
            fake_history(120.0, 0.15),
            fake_history(130.0, 0.25),
            fake_history(140.0, 0.35),
            fake_history(150.0, 0.45),
        ];

        // let errors = compare_growth(&measured, &history, false).unwrap();

        // let mut e = errors.into_iter();

        // assert!((e.next().unwrap() - 0.1054).abs() < 1e-3);
        // assert!((e.next().unwrap() - -0.09531).abs() < 1e-3);
    }

    fn fake_history(block: f64, a: f64) -> grow::History {
        let c = 0.2;
        let geometry = CrackGeometry {
            a: Some(CrackFront{
                length: a,
                angle: 0.0,
            }),
            c: Some(CrackFront{
                length: c,
                angle: FRAC_PI_2,
            }),
            ratio: Some(a / c),
        };

        grow::History {
            block,
            kmax: StressIntensity {a: Some(10.0), c: Some(20.0)},
            kmin: StressIntensity::default(),
            dk: StressIntensityDelta {a: Some(5.0), c: Some(10.0)},
            stress: 1.0,
            cycle: cycle::Cycle {
                min: tag::Tag::new(5.0, 0),
                max: tag::Tag::new(10.0, 0),
            },
            beta: BetaResult {a: Some(0.1), c: Some(0.2)},
            length_deltas: CrackLengthDelta {a: Some(0.1), c: Some(0.3)},
            geometry,
        }
    }

    #[test]
    fn normalise_returns_correct_scale_when_largest_at_start() {
        let largest = 3.0;
        let mut data = vec![largest, 1.0, 2.0];
        let scale = _normalise(&mut data);
        let difference = largest - scale;

        assert!(difference.abs() <= EPSILON);
    }

    #[test]
    fn normalise_returns_correct_scale_when_largest_at_end() {
        let largest = 3.0;
        let mut data = vec![1.0, 2.0, largest];
        let scale = _normalise(&mut data);
        let difference = largest - scale;

        assert!(difference.abs() <= EPSILON);
    }

    #[test]
    fn normalise_returns_correct_scale_when_largest_is_negative() {
        let largest = -3.0;
        let mut data = vec![largest, 1.0, 2.0];
        let scale = _normalise(&mut data);
        let difference = largest.abs() - scale;

        assert!(difference.abs() <= EPSILON);
    }

    //#[test]
    fn _test_range() {
        
        let mut options = options::get_default_options();
        options.seq_infile = "FALSTAFF.txt".to_string();
        options.seq_mods.reorder = true;
        options.scale = 200.0;
        options.beta = "qeft-newman84".to_string();
        options.a = 7e-5;
        options.c = 7e-5;
        options.a_limit = 0.06;
        options.c_limit = 0.06;
        options.component.forward = 0.025;
        options.component.sideways = 0.00635;
        options.dadn = "paris".to_string();
        options.params = BTreeMap::new();
        options.params.insert(ParameterLabel::c, 1e-9);
        options.params.insert(ParameterLabel::m, 3.0);
        options.optimise.file = "optimCoupons.txt".to_string();
        options::read_all_files(&mut options);

        // let mut fatigue_test = options::builder::grow::get_fatigue_test(&options);
        // let mut samples = options::builder::optimise::get_all_samples(&options);
        
        // let mut optimisation = options::builder::optimise::get_optimisation(&options);
        let mut samples = options::builder::optimise::get_all_samples(&options);
        // optimisation.run(&mut samples);
        
        
        

        // fatigue_test.run();
        // println!("{:#?}", fatigue_test.get_history().last().unwrap());

        let mut params = vec![1e-9, 2.5];

        for i in 0..=10 {
            params[1] = 2.95 + (i as f64 / 100.0);
            println!("{}, {}", params[1], sum_prediction_error(&mut samples, &params, false));
        }
    }
}
