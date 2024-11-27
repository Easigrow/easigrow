use self::rayon::prelude::*;
use crate::factors;
use crate::numbers::NonNan;
use crate::optimise;
use std::f64;

use super::Optimisation;

extern crate rayon;

pub struct Sweep {
    /// Parameters to be optimised
    parameters: Vec<f64>,
    /// Factors used to scale each parameter
    factors: Vec<f64>,
    /// Printable result
    result: String,
}

impl Sweep {
    pub fn new(parameters: Vec<f64>, factors: Vec<f64>) -> Self {
        let result = "N/A".to_string();
        Self {
            parameters,
            factors,
            result,
        }
    }

    /// Perform a brute force sweep over all permutations of the parameters.
    ///
    /// Calculates the error function for each permutation of the
    /// factors.  This function can be used as a brute force optimisation to search
    /// a space for the best fit.  Or, it could be used for a one off error
    /// function evaluation for all of the crack files that are to be matched.
    pub fn run(&mut self, samples: &mut [optimise::Sample], use_log_crack_depth: bool) {
        println!("Performing a brute force sweep. This may take some time...");
        let sweep_factors = factors::permutations(&self.factors, self.parameters.len());
        println!("Sweep: there are {} combinations", sweep_factors.len());
        println!("Sweep: parameters {:?}", &self.parameters);

        let mut data = Vec::new();
        for factors in &sweep_factors {
            data.push((factors, samples.to_owned()));
        }

        let mut results = Vec::with_capacity(sweep_factors.len());
        data.par_iter_mut()
            .map(|(normalised_factors, samples)| {
                let scaled_factors = normalised_factors
                    .iter()
                    .zip(self.parameters.iter())
                    .map(|(&f, &p)| f * p)
                    .collect::<Vec<_>>();
                let single_error =
                    optimise::prediction_error(&scaled_factors, samples, use_log_crack_depth);
                let total_single_error = single_error
                    .iter()
                    .map(|x| optimise::rms(x).powi(2))
                    .sum::<f64>()
                    .sqrt();
                (scaled_factors, total_single_error)
            })
            .collect_into_vec(&mut results);

        for (f, r) in &results {
            println!("{:?} {:?}", f, r);
        }

        // print out the results for the smallest error
        let (best_scaled_factors, smallest_error) = results
            .into_iter()
            .min_by_key(|(_factor, error)| NonNan::new(error.abs()))
            .unwrap();

        println!("Best result from sweep:");
        println!("    Total Error: {:?}", smallest_error);
        println!("    Scaled factors: {:?}", best_scaled_factors);

        // copy back the best parameters
        let n = self.parameters.len();
        self.parameters.clone_from_slice(&best_scaled_factors[..n]);

        self.result = samples[0].fatigue_test.printable_parameters(&self.parameters)
    }
}

impl Optimisation for Sweep {
    fn run(&mut self, samples: &mut [optimise::Sample], use_log_crack_depth: bool) {
        self.run(samples, use_log_crack_depth)
    }

    fn get_parameters(&self) -> &[f64] {
        &self.parameters
    }

    fn get_result(&self) -> String {
        self.result.to_owned()
    }
}
