use super::dadn::*;
use crate::options::{EasiOptions, GrowthMethod};
use fatigue::{
    beta,
    cycle::{cycles_from_sequence, process_cycle_mods, process_seq_mods},
    dadn::ParameterLabel,
    grow::{
        self,
        stress_constraint::{StressConstraints, VariableAlpha},
        CrackFront, CrackGeometry, CrackLength, FatigueTest, Grow,
    },
    tag::Tag,
};
use log::error;

/// Returns the fatigue test as well as any unclosed cycles
pub fn get_fatigue_test(options: &EasiOptions) -> (Box<dyn FatigueTest + Send + Sync>, Vec<Tag>) {
    if options.scale == 0.0 {
        error!(
            "Error: The sequence scale factor is 0. You need to set the scale factor
        (i.e. load or stress level) in order to perform a crack growth calculation.
        Try\n easigrow --help"
        );
        std::process::exit(1);
    }

    let dadn_params = match get_dadn_params(options) {
        Ok(result) => result.0,
        Err(why) => {
            error!("Error: {}", why);
            std::process::exit(1)
        }
    };
    let dadn = get_dadn(options, &dadn_params);

    let geometry = CrackGeometry {
        a: Some(CrackFront {
            length: options.a,
            angle: options.a_angle,
        }),
        c: Some(CrackFront {
            length: options.c,
            angle: options.c_angle,
        }),
        ratio: Some(options.a / options.c),
    };
    let beta = beta::get_beta_fn(
        &options.beta,
        &options.component,
        geometry,
        options.beta_doi.clone(),
    );

    let initial = CrackLength {
        a: Some(options.a),
        c: Some(options.c),
    };

    let limit = CrackLength {
        a: Some(options.a_limit),
        c: Some(options.c_limit),
    };

    let mut grow = Grow {
        cycles: vec![],
        component: options.component.clone(),
        scale: options.scale,
        limit,
        block_limit: options.block_limit,
        termination_tests: options.termination_tests.clone(),
        keep_every: options.optimise.keep_every,
    };

    match options.method {
        GrowthMethod::Easigrow => {
            let sequence = process_seq_mods(&options.sequence, &options.seq_mods);

            let (cycles, unclosed) = if options.cycle_infile.is_empty() {
                cycles_from_sequence(&sequence, &options.cycle_method)
            } else {
                (options.cycles.clone(), Vec::new())
            };

            // process all the modifications to the cycles
            grow.cycles = process_cycle_mods(&cycles, &options.cycle_mods);

            (Box::new(grow::easigrow::Easigrow::new(grow, initial, dadn, beta))
                as Box<dyn FatigueTest + Send + Sync>,
            unclosed)
        }

        GrowthMethod::Willenborg => {
            let sequence = process_seq_mods(&options.sequence, &options.seq_mods);
            let stress_constraints = get_stress_constraints(options);
            let deltak_th = get_deltak_th(options);

            let (cycles, unclosed) = if options.cycle_infile.is_empty() {
                cycles_from_sequence(&sequence, &options.cycle_method)
            } else {
                (options.cycles.clone(), vec![])
            };

            grow.cycles = process_cycle_mods(&cycles, &options.cycle_mods);

            (Box::new(grow::willenborg::Willenborg::new(
                grow,
                initial,
                dadn,
                beta,
                stress_constraints.alpha,
                stress_constraints.variable_alpha,
                deltak_th,
                options.r_so,
            )) as Box<dyn FatigueTest + Send + Sync>,
            unclosed)
        }
    }
}

pub fn get_deltak_th(options: &EasiOptions) -> f64 {
    // Get a deltak threshold value from the dadn parameters. If one isn't found,
    // as in the case of equations that don't use it, check for a value in the
    // separate CLI option.
    match options.params.get(&ParameterLabel::deltak_th) {
        Some(deltak_th) => *deltak_th,
        None => match options.deltak_th {
            Some(deltak_th) => deltak_th,
            None => {
                error!("This growth method requires a deltak_th value. Please specify one with the --deltak_th option.");
                std::process::exit(1)
            }
        },
    }
}

fn get_stress_constraints(options: &EasiOptions) -> StressConstraints {
    let mut result = StressConstraints {
        variable_alpha: None,
        alpha: None,
    };

    // If we've been given a variable alpha, use it
    if let Some(data) = &options.variable_alpha {
        let calculation_mode = match options.variable_alpha_mode {
            crate::options::VariableAlphaMode::GrowthRate => {
                grow::stress_constraint::CalculationMode::GrowthRate
            }
            crate::options::VariableAlphaMode::LogGrowthRate => {
                grow::stress_constraint::CalculationMode::LogGrowthRate
            }
            crate::options::VariableAlphaMode::CrackLength => {
                grow::stress_constraint::CalculationMode::CrackLength
            }
            crate::options::VariableAlphaMode::LogCrackLength => {
                grow::stress_constraint::CalculationMode::LogCrackLength
            }
        };

        match VariableAlpha::from_map(data, calculation_mode) {
            Ok(variable_alpha) => result.variable_alpha = Some(variable_alpha),
            Err(msg) => {
                error!("Error creating variable constraint: {}", msg);
                std::process::exit(1);
            }
        }
    } else {
        // Otherwise, use a constant value
        result.alpha = Some(options.alpha);
    };

    result
}
