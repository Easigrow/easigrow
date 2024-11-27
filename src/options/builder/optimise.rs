use std::{
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
    process,
};

use fatigue::grow::{willenborg::Willenborg, OptimisableParameter};
use log::error;

use crate::{
    optimise::{
        nelder::Nelder,
        particle_swarm::{Bound, ParticleSwarm},
        Optimisation, Sample, sweep::Sweep,
    },
    options::{self, clap::get_options_clap, EasiOptions},
};

use super::grow::get_fatigue_test;

/// Read all the option lines from the optimisation file.
pub fn get_all_options(main_options: &EasiOptions) -> Vec<EasiOptions> {
    let path = Path::new(&main_options.optimise.file);
    let handle = match File::open(path) {
        Err(why) => {
            error!("Error: could not open file '{}': {}.", path.display(), why);
            process::exit(1)
        }
        Ok(file) => file,
    };

    let optimfile = BufReader::new(handle);

    // create a vector of the process options for each optimisation
    let mut all_options: Vec<EasiOptions> = Vec::new();

    // read in each line of the optimisation file and store as option arguments
    for line in optimfile.lines() {
        let mut options = main_options.clone();

        match line {
            Ok(line) => {
                if main_options.verbosity == options::Verbosity::Verbose {
                    println!("matchfile: {:?}", line);
                }
                get_options_clap(&line, &mut options);
            }
            Err(e) => {
                error!(
                    "Error: problem in reading line from the optimisation file: {}.",
                    e
                );
                process::exit(1);
            }
        }

        if options.crack_infile.is_empty() {
            error!(
                "Error: The --crack_infile option is missing from the list of
            parameters in the optimisation file. The crackfile option
            is necessary for optimisation because it supplies the target
            crack growth curve."
            );
            process::exit(1);
        }

        if !options.cycle_infile.is_empty() && !options.seq_infile.is_empty() {
            error!("Error: you have specified a sequence file '{}' as well as a cycle file '{}'. Specify only one.",
                     options.seq_infile, options.cycle_infile);
            process::exit(2)
        }

        // make sure we are using the main dadn model in every sub model
        options.dadn.clone_from(&main_options.dadn);
        println!("crackfile: {}", options.crack_infile);
        options::read_all_files(&mut options);

        all_options.push(options);
    }
    all_options
}

pub fn get_all_samples(options: &EasiOptions) -> Vec<Sample> {
    let all_options = get_all_options(options);
    let mut samples = Vec::new();

    for option in all_options {
        let (fatigue_test, _unclosed) = get_fatigue_test(&option);
        let sample = Sample {
            fatigue_test,
            measurements: option.fracto.clone(),
            weight: option.crack_weight,
        };

        samples.push(sample);
    }

    samples
}

pub fn get_optimisation(options: &EasiOptions) -> Box<dyn Optimisation> {
    match options.optimise.method {
        options::OptimMethod::Sweep => Box::new(Sweep::new(
            get_parameters(options),
            options.optimise.sweep.to_owned(),
        )),
        options::OptimMethod::Nelder => Box::new(Nelder::new(
            options.optimise.nelder_params,
            get_parameters(options),
            options.optimise.tol,
            options.optimise.maxiter,
        )) as Box<dyn Optimisation>,
        options::OptimMethod::Particle => Box::new(ParticleSwarm::new(
            options.optimise.swarm_size,
            get_parameters(options),
            options.optimise.maxiter,
            options.optimise.min_fitness,
            get_bounds(options),
            options.optimise.bounds_constraints,
        )) as Box<dyn Optimisation>,
    }
}

fn get_parameters(options: &EasiOptions) -> Vec<f64> {
    let method_parameters = get_method_optimisable_parameter_values(options);
    let mut dadn_parameters = match super::dadn::get_dadn_params(options) {
        Ok(result) => result.0,
        Err(why) => {
            error!("Error: {}", why);
            std::process::exit(1)
        }
    };

    let mut result = vec![];

    for (_, value) in method_parameters {
        result.push(value);
    }

    result.append(&mut dadn_parameters);

    result
}

fn get_bounds(options: &EasiOptions) -> Vec<Bound> {
    // Add the method parameters first
    let mut bounds = get_method_optimisable_bounds(options);

    // Each parameter will have a specified min and max value. If either of these
    // values has not been provided, then the initial parameter value will be used.
    // This means if both min and max are omitted for a parameter, it will remain
    // fixed throughout the optimisation.
    for (parameter_label, value) in &options.params {
        let default = (*value, *value);
        let values = *options
            .optimise
            .range_dadn
            .get(parameter_label)
            .unwrap_or(&default);

        bounds.push(Bound::from_tuple(values));
    }

    bounds
}

fn get_method_optimisable_parameters(options: &EasiOptions) -> Vec<OptimisableParameter> {
    match options.method {
        options::GrowthMethod::Easigrow => vec![],
        options::GrowthMethod::Willenborg => Willenborg::optimisable_parameters(options.variable_alpha.is_some()),
    }
}

fn get_method_optimisable_parameter_values(
    options: &EasiOptions,
) -> Vec<(OptimisableParameter, f64)> {
    use OptimisableParameter::*;

    // Collect all the disparate values. Must maintain original order of values as
    // they have been given because the growth method depends on the order.
    let mut result: Vec<(OptimisableParameter, f64)> = Vec::new();
    for parameter in get_method_optimisable_parameters(options) {
        let value = match parameter {
            rso => options.r_so,
            alpha => options.alpha,
            deltak_th => super::grow::get_deltak_th(options),
            alpha_min => *options.variable_alpha.as_ref().unwrap().get(alpha_min.text()).unwrap(),
            alpha_max => *options.variable_alpha.as_ref().unwrap().get(alpha_max.text()).unwrap(),
            k => *options.variable_alpha.as_ref().unwrap().get(k.text()).unwrap(),
            x0 => *options.variable_alpha.as_ref().unwrap().get(x0.text()).unwrap(),
        };

        result.push((parameter, value));
    }

    result
}

fn get_method_optimisable_bounds(options: &EasiOptions) -> Vec<Bound> {
    let mut bounds: Vec<Bound> = Vec::new();

    let method_parameters = get_method_optimisable_parameter_values(options);

    for (optimisable_parameter, value) in method_parameters {
        let default = (value, value);
        let values = *options
            .optimise
            .range_method
            .get(&optimisable_parameter)
            .unwrap_or(&default);

        bounds.push(Bound::from_tuple(values));
    }

    bounds
}
