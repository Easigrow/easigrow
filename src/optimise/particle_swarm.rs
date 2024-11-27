//! Numerical optimisation using a Particle Swarm algorithm
//!
//! References:
//!     Maurice Clerc. Standard Particle Swarm Optimisation. 2012. hal-00764996
//!         https://hal.archives-ouvertes.fr/hal-00764996
//!     Freitas, D., Lopes, L. and Morgado-Dias, F., 2020. Particle Swarm Optimisation: A Historical Review Up to the Current Developments.
//!         https://www.mdpi.com/1099-4300/22/3/362

// Overview of the general form of the algorithm
// Initialise the swarm
//      Pick a random position
//      Pick a random velocity
//      Compute the fitness, setting previous best to this initial position/fitness
//      Create the neighbourhoods
// Iterate
//      For each particle:
//          Compute the new velocity
//          Move the particle using the newly calculated velocity
//          *Optional: Apply a confinement method
//          If new position is better, update previous best position and fitness
//      Stop if:
//          Error is within acceptable range, or
//          Maximum iterations has been reached

use rand::{distributions::Uniform, prelude::*};
use rayon::prelude::*;
use std::sync::{Arc, Mutex};
use std::time::SystemTime;

use super::{sum_prediction_error, Optimisation, Sample};

#[derive(Debug, Clone)]
pub struct Bound {
    min: f64,
    max: f64,
}

#[derive(Clone, Debug)]
struct Evaluation {
    position: Vec<f64>,
    fitness: f64,
}

#[derive(Debug)]
struct Particle {
    current: Evaluation,
    best: Evaluation,
    velocity: Vec<f64>,
    neighbours: Vec<usize>,
    neighbourhood_best: Evaluation,
    is_neighbourhood_best: bool,
    // The below are growth specific data. They are currently stored here to
    // facilitate concurrent evaluations. If a single samples reference were 
    // shared between all particles, it would need to be behind a mutex. Doing
    // this would mean evaluations are no longer concurrent - making an 
    // optimisation take far too long. It would be good if they could be carved 
    // off so the implementation was more generic, though this may not be possible.
    samples: Vec<Sample>,
    use_log_crack_depth: bool,
}

#[derive(Debug)]
pub struct ParticleSwarm {
    size: usize,
    particles: Vec<Particle>,
    best: Evaluation,
    max_iterations: usize,
    min_fitness: f64,
    bounds: Vec<Bound>,
    bounds_constraints: bool,
    start_time: SystemTime,
}

impl Bound {
    /// Returns None if min > max
    pub fn _new(min: f64, max: f64) -> Option<Self> {
        if min > max {
            None
        } else {
            Some(Self {
                min,
                max,
            })
        }
    }

    pub fn from_tuple(values: (f64, f64)) -> Self {
        let min;
        let max;

        if values.0 < values.1 {
            min = values.0;
            max = values.1;
        } else {
            min = values.1;
            max = values.0;
        }

        Self {
            min,
            max,
        }
    }
}

impl Particle {
    fn new(position: Vec<f64>, velocity: Vec<f64>, samples: Vec<Sample>, use_log_crack_depth: bool,) -> Particle {
        Particle {
            current: Evaluation {
                position: position.clone(),
                fitness: f64::INFINITY,
            },
            velocity,
            best: Evaluation {
                position: position.clone(),
                fitness: f64::INFINITY,
            },
            neighbours: vec![],
            neighbourhood_best: Evaluation {
                position: position.clone(),
                fitness: f64::INFINITY,
            },
            is_neighbourhood_best: false,
            samples,
            use_log_crack_depth,
        }
    }

    fn evaluate(&mut self) -> f64 {
        // let mut samples = &mut *self.samples.lock().unwrap(); 
        let result = sum_prediction_error(&mut self.samples, self.current.position.as_slice(), self.use_log_crack_depth);
        self.current.fitness = result;
        result
    }
}

impl ParticleSwarm {
    // Denoted 'K' in the adaptive random topology
    const NEIGHBOURS: u32 = 3;

    pub fn new(
        size: usize,
        parameters: Vec<f64>,
        max_iterations: usize,
        min_fitness: f64,
        bounds: Vec<Bound>,
        bounds_constraints: bool,
    ) -> ParticleSwarm {
        ParticleSwarm {
            size,
            particles: vec![],
            best: Evaluation {
                position: parameters,
                fitness: f64::INFINITY,
            },
            max_iterations,
            min_fitness,
            bounds,
            bounds_constraints,
            start_time: SystemTime::now(),
        }
    }

    fn initialise(&mut self, samples: Vec<Sample>, use_log_crack_depth: bool) {
        let mut particles = Vec::new();
        self.start_time = SystemTime::now();

        for i in 0..self.size {
            let mut rng = rand::thread_rng();
            // Pick a random position and velocity
            let mut position: Vec<f64> = Vec::new();
            let mut velocity: Vec<f64> = Vec::new();

            if i == 0 {
                // Create an "exact" particle using the initial parameter values
                position.clone_from(&self.best.position);

                for (i, bound) in self.bounds.iter().enumerate() {
                    let velocity_dimension =
                        rng.gen_range((bound.min - position[i])..=(bound.max - position[i]));
                    velocity.push(velocity_dimension);
                }
            } else {
                for bound in &self.bounds {
                    let position_dimension = rng.gen_range(bound.min..=bound.max);
                    let velocity_dimension = rng.gen_range(
                        (bound.min - position_dimension)..=(bound.max - position_dimension),
                    );
                    position.push(position_dimension);
                    velocity.push(velocity_dimension);
                }
            }

            // Create the particle
            let particle = Particle::new(position, velocity, samples.clone(), use_log_crack_depth);

            particles.push(particle);
        }

        // Do an initial evaluation on all the particles
        particles.par_iter_mut().for_each(|particle| {
            particle.best.fitness = particle.evaluate();
        });

        self.particles = particles;
        self.set_all_neighbourhoods();
    }

    pub fn run(&mut self) -> f64 {
        // TODO: Fix this cloning. It's being done because the loop over
        // the particles is complaining about borrowing twice
        let bounds_constraints = self.bounds_constraints;
        let bounds = self.bounds.clone();
        let min_fitness = self.min_fitness;

        let mut best_iteration = 0;
    
        // Initialise parameter values
        let w = 1.0 / (2.0 * 2.0_f64.ln());
        let c = 0.5 + 2.0_f64.ln();
        let min_fitness_reached = Arc::new(Mutex::new(false));
    
        for i in 1..=self.max_iterations {
            if *min_fitness_reached.lock().unwrap() {
                break;
            }
    
            self.print_iteration_time(i);
    
            // Iterate over each particle, computing new positions and fitness values
            self.particles.par_iter_mut().for_each(|particle| {
                let mut rng = rand::thread_rng();
    
                // Update velocity and position
                #[allow(clippy::needless_range_loop)]
                for dimension in 0..particle.velocity.len() {
                    let mut result = w * particle.velocity[dimension];
                    result += rng.gen_range(0.0..=c)
                        * (particle.best.position[dimension] - particle.current.position[dimension]);
                    if !particle.is_neighbourhood_best {
                        result += rng.gen_range(0.0..=c)
                            * (particle.neighbourhood_best.position[dimension]
                                - particle.current.position[dimension]);
                    }
                    particle.velocity[dimension] = result;
                    particle.current.position[dimension] += result;
    
                    // Confinement
                    if bounds_constraints {
                        if particle.current.position[dimension] < bounds[dimension].min {
                            particle.current.position[dimension] = bounds[dimension].min;
                            particle.velocity[dimension] *= -0.5;
                        } else if particle.current.position[dimension] > bounds[dimension].max {
                            particle.current.position[dimension] = bounds[dimension].max;
                            particle.velocity[dimension] *= -0.5;
                        }
                    }
                }
    
                // Evaluate new position
                particle.evaluate();
                if particle.current.fitness < particle.best.fitness {
                    particle.best.fitness = particle.current.fitness;
                    particle.best.position.clone_from(&particle.current.position);
                }

                if particle.best.fitness < min_fitness {
                    *min_fitness_reached.lock().unwrap() = true;
                    println!(
                        "Minimum fitness reached: {}",
                        particle.best.fitness
                    );
                }
            });
    
            self.update_neighbourhood_best();
            if self.update_global_best() {
                best_iteration = i;
                println!("New global best: {:?}", self.best);
            }
        }
    
        println!("\nOptimisation complete");
        // println!("Optimised parameters: {:?}", relabel_parameters(&self.best.position, &main_options.dadn).unwrap());
        println!("Best fitness: {}", self.best.fitness);
        println!("Achieved on iteration: {}", best_iteration);
        self.print_final_time();

        self.best.fitness
    }

    fn set_all_neighbourhoods(&mut self) {
        for i in 0..self.particles.len() {
            self.set_neighbourhood(i);
        }
    }

    fn set_neighbourhood(&mut self, particle_index: usize) {
        let mut neighbours: Vec<usize> = Vec::new();

        // A particle always informs at least itself
        neighbours.push(particle_index);

        // Pick randomly from the swarm. This method will mean any particle will be
        // informed by at least itself and at most NEIGHBOURS + 1 particles (including itself).
        let mut rng = rand::thread_rng();
        let range = Uniform::from(0..self.size);
        for _ in 0..ParticleSwarm::NEIGHBOURS {
            neighbours.push(range.sample(&mut rng));
        }

        self.particles[particle_index].neighbours = neighbours;
    }

    fn update_neighbourhood_best(&mut self) {
        // TODO: Investigate if there is a more idiomatic way to achieve this. At the moment
        // it looks very C-style instead of Rust.

        for i in 0..self.particles.len() {
            let mut is_neighbourhood_best = true;
            let mut best = &self.particles[i].best;

            for neighbour_index in &self.particles[i].neighbours {
                if self.particles[*neighbour_index].best.fitness < best.fitness {
                    best = &self.particles[*neighbour_index].best;
                    is_neighbourhood_best = false;
                }
            }

            self.particles[i].neighbourhood_best = best.clone();
            self.particles[i].is_neighbourhood_best = is_neighbourhood_best;
        }
    }

    fn update_global_best(&mut self) -> bool {
        let mut has_updated = false;

        for particle in &self.particles {
            if particle.best.fitness < self.best.fitness {
                self.best = particle.best.clone();
                has_updated = true;
            }
        }

        // Adaptive random topology re-allocates neighbourhoods when the
        // global best has not improved after an iteration
        if !has_updated {
            for i in 0..self.particles.len() {
                self.set_neighbourhood(i);
            }

            self.update_neighbourhood_best();
        }

        has_updated
    }

    fn print_iteration_time(&self, iteration: usize) {
        println!("\nIteration: {} (of {})", iteration, self.max_iterations);
        match self.start_time.elapsed() {
            Ok(elapsed) => {
                let elapsed_millis = elapsed.as_millis();
                let elapsed = elapsed.as_secs();
                let secs = elapsed % 60;
                let mins = (elapsed / 60) % 60;
                let hours = elapsed / 3600;
                println!("Time elapsed : {}h {:02}m {:02}s", hours, mins, secs);

                let per_iteration_estimate = elapsed_millis / iteration as u128;
                let remaining = (per_iteration_estimate * (self.max_iterations + 1 - iteration) as u128) / 1000;
                let secs = remaining % 60;
                let mins = (remaining / 60) % 60;
                let hours = remaining / 3600;
                println!(
                    "Time remaining (approx.): {}h {:02}m {:02}s",
                    hours, mins, secs
                );
            }
            Err(e) => {
                println!("Timing error: {:?}", e);
            }
        }
    }

    fn print_final_time(&self) {
        match self.start_time.elapsed() {
            Ok(elapsed) => {
                let elapsed = elapsed.as_secs();
                let secs = elapsed % 60;
                let mins = (elapsed / 60) % 60;
                let hours = elapsed / 3600;
                println!("Total time elapsed: {}h {:02}m {:02}s", hours, mins, secs);
            }
            Err(e) => {
                println!("Error: {:?}", e);
            }
        }
    }
}

impl Optimisation for ParticleSwarm {
    fn get_parameters(&self) -> &[f64] {
        &self.best.position
    }

    fn run(&mut self, samples: &mut [Sample], use_log_crack_depth: bool) {
        self.initialise(samples.to_vec(), use_log_crack_depth);
        self.run();
    }

    fn get_result(&self) -> String {
        // We assume that there is at least 1 particle and at least 1 sample
        // We also assume that all samples use the same growth method
        self.particles[0].samples[0].fatigue_test.printable_parameters(&self.best.position)
    }
}
