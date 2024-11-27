//! Numerical optimisation using the Nelder-Mead algorithim
//!
//! Translation of the pure Python/Numpy implementation of the Nelder-Mead algorithm.
//! Reference: <https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method>

use std::process;

// We use this crate to allow us to add and multiply vectors to keep it in line with python.
// similar to nalgebra but nalgrebra requires way too many subcrates
use crate::vector::MVector;
use log::{warn, error};

use super::{sum_prediction_error_factorised, Optimisation, Sample};

/// A point on the simplex.
#[derive(Debug, Clone)]
struct Eval {
    x: MVector<f64>,
    score: f64,
}

// finds the centroids of a list of points
fn centroid(simplex: &[Eval]) -> MVector<f64> {
    let n = simplex[0].x.len();
    let zero = vec![0.0_f64; n];

    //    let x0 = MVector::from_slice(n, &zero);
    let x0 = MVector::from_row_slice(n, &zero);
    simplex.iter().fold(x0, |sum, r| sum + r.x.clone()) / simplex.len() as f64
}

// Keep track of each type of modification.
#[derive(Debug)]
enum Operation {
    Reflection,
    Expansion,
    Contraction,
    Reduction,
}

/// Nelder-Mead specific optimisation parameters.
#[derive(Debug, Clone, Copy)]
pub struct Parameters {
    /// step increment in each direction from the starting point
    pub step: f64,
    /// reflection factor (relects the worst point through the centroid)
    pub alpha: f64,
    /// expansion factor (extends the reflected point)
    pub gamma: f64,
    /// contraction factor (away from the worst point)
    pub rho: f64,
    /// shrinking factor (around the best point)
    pub sigma: f64,
}

#[derive(Debug)]
pub struct Nelder {
    /// Nelder-Mead specific parameters
    nelder_parameters: Parameters,
    /// Parameters to be optimised
    parameters: Vec<f64>,
    //TODO: This should probably go into nelder_parameters
    converge_tol: f64,
    /// Maximum number of iterations to run
    max_iter: usize,
    /// Printable result
    result: String,
}

impl Parameters {
    pub fn new(values: &[f64]) -> Parameters {
        Parameters {
            step: values[0],
            alpha: values[1],
            gamma: values[2],
            rho: values[3],
            sigma: values[4],
        }
    }

    pub fn default() -> Parameters {
        Parameters {
            step: 0.1,
            alpha: 1.0,
            gamma: 2.0,
            rho: 0.5,
            sigma: 0.5,
        }
    }
}

impl Nelder {
    pub fn new(
        nelder_parameters: Parameters,
        parameters: Vec<f64>,
        converge_tol: f64,
        max_iter: usize,
    ) -> Self {
        let result = "N/A".to_string();
        Self {
            nelder_parameters,
            parameters,
            converge_tol,
            max_iter,
            result,
        }
    }

    /// Nelder-Mead non-linear optimisation
    ///
    /// This routine works best if each of the parameters being optimised are
    /// roughly the same size.  If this is not the case then they should
    /// be normalised to ensure they are.
    fn run<F>(&mut self, mut objective_function: F) -> f64
    where F: FnMut(&[f64]) -> f64 {
        check_nelder_limits(&self.nelder_parameters);
        // let factors = vec![1.0; self.parameters.len()];

        let x_start = MVector::from_row_slice(self.parameters.len(), &self.parameters);
        let mut results = Vec::new();
        let mut ops = Vec::new();

        // results.push(Result {x: x_start.clone(), score: f(x_start.as_ref())});
        results.push(Eval {
            x: x_start.clone(),
            score: objective_function(x_start.as_slice()),
        });

        // iniitalise the simplex by taking a step in each direction.
        for i in 0..x_start.len() {
            let mut x_init = x_start.clone();
            // just have to peer inside vec here cause I can't make it work otherwise
            // x_init.m[i] *= 1.0 + self.nelder_parameters.step;
            // if x_init.m[i] == 0.0 {
            //     x_init.m[i] = self.nelder_parameters.step;
            // }
            x_init.m[i] += self.nelder_parameters.step;
            // results.push(Eval {x: x_init.clone(), score: f(x_init.as_ref())});
            results.push(Eval {
                x: x_init.clone(),
                score: objective_function(x_init.as_slice()),
            });
        }
        println!("Nelder: starting result {:?}", results);

        let mut prev_best = results[0].score;
        let mut iter = 0;
        let n = results.len();

        loop {
            // Check if exceeding the iteration limit.
            if iter >= self.max_iter {
                warn!("***Warning: The optimisation has failed to converge within the specified maximum iteration limit {}. 
    The answer may not be optimum. Try increasing the limit.", self.max_iter);
                break;
            }
            iter += 1;

            results.sort_by(|a, b| a.score.partial_cmp(&b.score).unwrap());
            println!("Nelder: {} {}. best: {}", iter, n, results[0].score);

            // check for convergence
            let change_tol = (results[0].score - results[n - 1].score).abs();
            if change_tol < self.converge_tol {
                println!("Nelder: Success. Converged tol {}", change_tol);
                break;
            } else {
                println!("Nelder: convergence tolerance not reached {}", change_tol);
                println!(
                    "Nelder: best {}, worst {}, prev_best {}",
                    results[0].score,
                    results[results.len() - 1].score,
                    prev_best
                );
            }

            prev_best = results[0].score;

            // calculate centro
            let x0 = centroid(&results[0..n - 1]);

            // reflection
            // if the reflected point is better than the second worst,
            // but not better than the best, then replace the worst point
            // with the reflected point
            let xr = x0.clone() + self.nelder_parameters.alpha * (x0.clone() - results[n - 1].x.clone());
            // let rscore = f(xr.as_ref());
            let rscore = objective_function(xr.as_slice());
            println!("Nelder: rscore: {}", rscore);

            if rscore >= results[0].score && results[n - 2].score > rscore {
                ops.push(Operation::Reflection);
                println!("Nelder: Including reflected point");

                results.pop().unwrap();
                results.push(Eval {
                    x: xr,
                    score: rscore,
                });
                continue;
            }

            // expansion
            // If the reflected point is the best point so far
            // then compute the expanded point
            if rscore < results[0].score {
                let xe = x0.clone() + self.nelder_parameters.gamma * (xr.clone() - x0.clone());
                // let escore = f(xe.as_ref());
                let escore = objective_function(xe.as_slice());
                println!("Nelder: expansion score: {}", escore);

                results.pop().unwrap();

                if escore < rscore {
                    ops.push(Operation::Expansion);
                    results.push(Eval {
                        x: xe,
                        score: escore,
                    });
                    continue;
                } else {
                    results.push(Eval {
                        x: xr,
                        score: rscore,
                    });
                    continue;
                }
            }

            // contraction
            // If the contracted point is better than the worst point,
            // replace the worst point with the contracted point
            let xc = x0.clone() + self.nelder_parameters.rho * (results[n - 1].x.clone() - x0);
            // let cscore = f(xc.as_ref());
            let cscore = objective_function(xc.as_slice());

            if cscore < results[n - 1].score {
                println!("contracting: {}", cscore);

                ops.push(Operation::Contraction);
                results.pop().unwrap();
                results.push(Eval {
                    x: xc,
                    score: cscore,
                });
                continue;
            }

            // reduction
            // For all but the best point, replace the point with
            // xi = x1 + sigma(xi - x)
            // This is a shrinking of the simplex around the best point
            println!("Nelder: reducing");

            ops.push(Operation::Reduction);
            for r in 1..results.len() {
                results[r].x =
                    results[0].x.clone() - self.nelder_parameters.sigma * (results[r].x.clone() - results[0].x.clone());
                // results[r].score = f(results[r].x.as_ref());
                results[r].score = objective_function(results[r].x.as_slice());
            }
        }

        println!("Nelder: Iterations: {}", iter);
        let count = count_nelder_ops(&ops);
        println!(
            "Nelder: Operations: reflection {}, expansion {}, contraction {}, reduction {}",
            count.reflection, count.expansion, count.contraction, count.reduction
        );
        for (i, &a) in results[0].x.clone().as_slice().to_vec().iter().enumerate() {
            self.parameters[i] = a;
        }
        println!("Nelder: Score: {}", results[0].score);

        results[0].score
    }
}

impl Optimisation for Nelder {
    fn get_parameters(&self) -> &[f64] {
        &self.parameters
    }

    fn get_result(&self) -> String {
        self.result.to_owned()
    }

    fn run(&mut self, samples: &mut [Sample], use_log_crack_depth: bool) {
        // number of variables to be optimised
        let n = self.parameters.len(); 

        // Calculate the total number of points to be matched across all
        // match files. This may be excessive and could be reduced to one
        // value for each measured crack curve.
        let m: usize = samples.iter().fold(0, |sum, optimisation| {
            let nrun = 1.0f64;
            sum + (optimisation.measurements.len() as f64 - nrun) as usize
        });

        println!("match_crack: n (variables) {} m (targets) {}", n, m);    

        if m < n {
            error!(
                "Error: Insufficient number of match points {} to optimise {} variables",
                m, n
            );
            process::exit(1);
        }

        // let scale = normalise(&mut self.parameters);
        let cloned_params = self.parameters.to_owned();
        self.parameters = vec![1.0; self.parameters.len()];
        
        let objective_function = |factors: &[f64]| sum_prediction_error_factorised(samples, &cloned_params, factors, use_log_crack_depth);

        self.run(objective_function);

        //rescale
        for (i, param) in self.parameters.iter_mut().enumerate() {
            *param *= cloned_params[i];
        }

        self.result = samples[0].fatigue_test.printable_parameters(&self.parameters);
    }
}

struct NelderSum {
    reflection: usize,
    expansion: usize,
    contraction: usize,
    reduction: usize,
}

// Count the number of each type of operation in the nelder search.
fn count_nelder_ops(ops: &[Operation]) -> NelderSum {
    ops.iter().fold(
        NelderSum {
            reflection: 0,
            expansion: 0,
            contraction: 0,
            reduction: 0,
        },
        |mut count, op| {
            match *op {
                Operation::Reflection => count.reflection += 1,
                Operation::Expansion => count.expansion += 1,
                Operation::Contraction => count.contraction += 1,
                Operation::Reduction => count.reduction += 1,
            };

            count
        },
    )
}

/// Check that the Nelder parameters are acceptable.
fn check_nelder_limits(params: &Parameters) {
    // recommended parameter limits (wikipedia)
    if params.gamma < 0.0 {
        warn!("***Warning: Wikipedia recommends using a Nelder-Mead value for gamma > 0.0, using {} .", params.gamma);
    }
    if { params.rho < 0.0 } | { params.rho > 0.5 } {
        warn!(
            "***Warning: Wikipedia recommends using a Nelder-Mead value 0.0 < rho < 0.5, using {}.",
            params.rho
        );
    }
    if params.sigma < 0.0 {
        warn!(
            "***Warning: Wikipedia recommends using a Nelder-Mead value for sigma > 0.0, using {}.",
            params.sigma
        );
    }
}



#[cfg(test)]
mod tests {
    use crate::optimise::Optimisation;
    use crate::optimise::nelder::Parameters;
    use crate::vector::MVector;
    use log::info;

    use super::centroid;
    use super::{Nelder, Eval};

    #[test]
    fn test_centroid() {
        let x = vec![
            Eval {
                x: MVector::from_row_slice(3, &[0.0, 0.0, 0.0]),
                score: 0.0,
            },
            Eval {
                x: MVector::from_row_slice(3, &[1.0, 2.5, 1.0]),
                score: 0.0,
            },
            Eval {
                x: MVector::from_row_slice(3, &[2.0, 3.5, 2.0]),
                score: 0.0,
            },
        ];

        let cent = centroid(&x);
        let ans = vec![1.0, 2.0, 1.0];

        for i in 0..3 {
            assert!((ans[i] - cent.m[i]).abs() < std::f64::EPSILON);
        }
    }

    #[test]
    fn test_nelder_sin() {
        let f = |x: &[f64]| (x[0].sin() * x[1].cos()) * (1.0 / (x[2].abs() + 1.0));

        let x = vec![
            Eval {
                x: MVector::from_row_slice(3, &[0.0f64, 0.0, 0.0]),
                score: 0.0,
            },
            Eval {
                x: MVector::from_row_slice(3, &[1.0f64, 1.5, 1.0]),
                score: 0.0,
            },
            Eval {
                x: MVector::from_row_slice(3, &[2.0f64, 2.5, 2.0]),
                score: 0.0,
            },
        ];
        info!("centroid: {:?}", centroid(&x));

        let x = vec![1.0f64, 2.0, 3.0];

        let f_x = f(&x);

        let nelder_parameters = Parameters {
            step: 0.1,
            alpha: 1.0,
            gamma: 2.0,
            rho: 0.5,
            sigma: 0.5,
        };

        let mut nelder = Nelder::new(
            nelder_parameters,
            x,
            1e-6,
            100,
        );

        let result = nelder.run(f);

        info!("Start result: {}", f_x);
        //info!("Optimimum result: {:?} at x: {:?}", result, x);
        assert!((result + 1.0).abs() < 0.01);
    }

    
    #[test]
    fn test_nelder_rosenbrock() {
        let x = vec![2.0f64, 3.0];

        let f_x = rosenbrock2d(&x);

        let nelder_parameters = Parameters {
            step: 0.1,
            alpha: 1.0,
            gamma: 2.0,
            rho: 0.5,
            sigma: 0.5,
        };

        let mut nelder = Nelder::new(
            nelder_parameters,
            x,
            1e-6,
            100,
        );

        let result = nelder.run(rosenbrock2d);

        let optimised_params = nelder.get_parameters();

        println!("Start result: {}", f_x);
        println!("Optimimum result: {:?} at x: {:?}", result, optimised_params);

        
        assert!((result).abs() < 0.01);
        assert!(optimised_params.iter().fold(0.0, |s, x| s + (x - 1.0).powi(2)) < 0.001);
    }

    fn rosenbrock2d(x: &[f64]) -> f64 {
        (1.0 - x[0]).powi(2) + 100.0 * (x[1] - x[0].powi(2)).powi(2)
    }

    #[test]
    fn test_nelder_himmelblau() {
        let x_all = vec![[5.0, -5.0], [-5.0, 5.0], [-5.0, -5.0], [0.0, 0.0]];
        let ans_all = vec![
            [3.584_428, -1.848_126],
            [-2.805_118, 3.131_312],
            [-3.779_310, -3.283_186],
            [3.0, 2.0],
        ];

        for (x, ans) in x_all.iter().zip(ans_all) {
            let f_x = himmelblau(x);
            let y = *x;

            let nelder_parameters = Parameters {
                step: 0.1,
                alpha: 1.0,
                gamma: 2.0,
                rho: 0.5,
                sigma: 0.5,
            };
    
            let mut nelder = Nelder::new(
                nelder_parameters,
                y.to_vec(),
                1e-10,
                100,
            );
    
            let result = nelder.run(himmelblau);

            let optimised_params = nelder.get_parameters();

            println!("Start result: {}", f_x);
            println!("Optimimum result: {:?} at x: {:?}", result, optimised_params);

            //assert!((result).abs() < 1e-5);
            let vec_error = optimised_params.iter()
                .zip(ans.iter())
                .fold(0.0, |s, (x, a)| s + (x - a).powi(2));
            println!(
                "vecs error {}, answer {:?} obtained {:?}",
                vec_error, ans, optimised_params
            );
            assert!(vec_error < 1e-5);
        }
    }

    fn himmelblau(x: &[f64]) -> f64 {
        // this has one local maximum at f(-0.270845, -0.923039) = 181.617
        // and four identical local minima at
        // f(3.0, 2.0) = 0.0
        // f(-2.805118, 3.131312) = 0.0
        // f(-3.779310, -3.283186) = 0.0
        // f(3.584428, -1.848126) = 0.0

        (x[0].powi(2) + x[1] - 11.0).powi(2) + (x[0] + x[1].powi(2) - 7.0).powi(2)
    }
}
