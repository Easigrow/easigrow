/// easiGrow
///
/// by Paul White (Nov 2014--2017)
/// written in rust (www.rust-lang.org)
///
/// A program to match crack growth predictions to measurements.
///
/// The program calculates fatigue crack growth rates and finds the
/// optimum parameters of a crack growth model to match predictions
/// with measurements.
///
/// **easiGrow** is a standalone program but most of the calculations
/// are done through calls to the associated **fatigue** library which
/// is included. The main program is for doing anything that
/// explicitly uses the command line flags inlcuding the optimisation
/// module. These flages are used to build the **EasiOptions** data
/// structure which is then used to generate the crack growth
/// history. The optimisation generates a crack growth curve which it
/// compares with a fractography file. It finds the error between
/// these measurements and tries to minimise the sum errors through
/// minimisation routines.
///
/// Currently, none of the models has a memory effect, so it is ok to
/// just start growing the crack from an iniital crack size that is
/// smaller than the initial fracto data. The struct `grow::CrackState`
/// also contains parameters that are passed along with the applied
/// loading _kmin_ and _kmax_, so any memory variables should be added to
/// this struct and will be availabe to be used by the _da/dn_ equation.
/// The simplest memory effect that is included in the `CrackState`
/// data is the plastic zone size, but there are no dadn equations
/// currently using this. The memory effect does not appear to be
/// strong in AA7050 material.
///
/// Think of the program flow as
///
/// 1. Read in data
/// 2. Filter the sequence (turning point, rainflow, risefall, deadband etc.) and convert to cycles
/// 3. Filter the list of cycles
/// 4. If required, optimise any parameters
/// 5. Perform a crack growth calculation
/// 6. Write out requested output

#[macro_use]
extern crate clap;
extern crate fatigue;
extern crate log;
extern crate env_logger;

use std::collections::HashMap;
use fatigue::dadn::relabel_parameters;
use fatigue::grow::{CrackFront, CrackGeometry};
use fatigue::{beta, cycle, dadn, fracto, grow, io};
use options::{clap::get_options_clap, builder, TerminatingOutput};
use fatigue::COMMENT;
use std::fs::File;
use std::path::Path;
use log::error;
use std::io::Write;

mod list;
mod optimise;
mod factors;
mod options;
mod numbers;
mod vector;

fn main() {
    env_logger::init();
    
    let mut options = options::get_default_options();
    get_options_clap("", &mut options);
    println!("{}easiGrow: version {}", COMMENT, crate_version!());
    println!("{}", COMMENT);
    if options.verbosity == options::Verbosity::Verbose {
        println!("{}Options: ", COMMENT);
        println!("{}", options);
    }
    
    if !options.cycle_infile.is_empty() && !options.seq_infile.is_empty() {
        error!("Error: you have specified a sequence file '{}' as well as a cycle file '{}'. Specify only one.",
                 options.seq_infile, options.cycle_infile);
        std::process::exit(2)
    }

    options::read_all_files(&mut options);

    // Check if we're just running to output data
    if options.output != TerminatingOutput::None {
        print_and_exit(&mut options);
    }
    
    let params = match builder::dadn::get_dadn_params(&options) {
        Ok(result) => {
            println!("{}", result.1);
            result.0
        }
        Err(why) => {
            error!("Error: {}", why);
            std::process::exit(1)
        }
    };
    if options.params.is_empty() && !options.is_dadn_from_file() {
        options.params = match relabel_parameters(&params, &options.dadn) {
            Ok(result) => result,
            Err(why) => {
                error!("{}", why);
                std::process::exit(1)
            },
        };
    }

    // Perform logic checks on the given inputs to ensure valid data
    // has been provided
    if let Err(why) = is_input_data_consistent(&options) {
        error!("{}", why);
        std::process::exit(1)
    }

    // Optimise the parameters to match the predicted crack growth
    // rates with the associated measured crack growth rates.
    if !options.optimise.file.is_empty() {
        let mut optimisation = builder::optimise::get_optimisation(&options);
        
        println!(
            "{}Now starting the optimisation with params {:?} ...",
            COMMENT, params
        );

        let mut samples = builder::optimise::get_all_samples(&options);
        optimisation.run(&mut samples, options.optimise.use_log_crack_depth);

        println!("{}...finished the optimisation. ", COMMENT);
        println!("{}The optimised parameters are: \n{}", COMMENT, optimisation.get_result());
        
        std::process::exit(0);
    }

    // Grow the crack
    let history_all = generate_crack_history(&options, &params);

    // Lastly, now that we've grown the crack, check if we need to
    // generate and write out a pseudo image.
    if !options.image.file.is_empty() {
        println!("Making a pseudo image...");
        if options.image.file.ends_with(".svg") {
            fracto::write_svg_pseudo_image(&history_all, &options.image);
            println!("Image written to file '{}'", options.image.file);
        } else {
            error!("Error: Currently easigrow can only generate svg. Please use a '.svg' suffix");
        }
    }
}

// Finally grow the crack with the current parameters which may have been optimised.

// We exit here if the scale has not been set. Otherwise we
// would go through and do a default calculation which confuses
// people if they just want to start the program to see how to get
// help.
fn generate_crack_history(options: &options::EasiOptions, params: &[f64]) -> Vec<grow::History> {
    let dadn_eqn = builder::dadn::get_dadn(options, params);
    println!("{}da/dN equation: {}", COMMENT, dadn_eqn);

    let (mut fatigue_test, _unclosed) = builder::grow::get_fatigue_test(options);

    // make a hash set of the lines that are required for output
    // let mut output_lines: BTreeSet<usize> = options.output_lines.iter().cloned().collect();

    // if there are no lines in the output then put in the line for the first cycle
//     if fatigue_test
//         .get_cycles()
//         .iter()
//         .filter(|c| output_lines.contains(&c.max.index) || output_lines.contains(&c.min.index))
//         .count() == 0
//     {
//         println!("output_lines {:?}", output_lines);
//         println!(
//             "
// Warning: There are no sequence lines in the cycle list and so there
//          will be no crack growth output. Consider closing up cycles
//          with re-order to use all sequence lines or include specific
//          sequence lines that are in the cycle. Meanwhile, the output will
//          be for the sequence line in the first cycle at line {}.",
//             options.cycles[0].max.index
//         );
//         output_lines.insert(options.cycles[0].max.index);
//     }

    let mut cycle_indexes = HashMap::new();
    for (i, cycle) in fatigue_test.get_cycles().iter().enumerate() {
        cycle_indexes.insert(cycle.max.index, i);
    }

    let (history, messages) = fatigue_test.run();

    if options.output_every.is_some() && !options.output_cycles.is_empty() {
        println!("Both --output_every and --output_cycles have been specified, only using --output_cycles")
    }

    grow::display_history_header(&options.output_vars);

    // Always display the first moment
    grow::display_history_line(history.first().unwrap(), &options.output_vars, &options.component, &cycle_indexes);
    
    if !options.output_cycles.is_empty() {
        let mut output_cycles = options.output_cycles.to_owned();
        output_cycles.sort();
        
        for moment in history.iter().skip(1) {
            if output_cycles.binary_search(cycle_indexes.get(&moment.cycle.max.index).unwrap()).is_ok() {
                grow::display_history_line(moment, &options.output_vars, &options.component, &cycle_indexes)
            }
        }
    } else {
        let output_every = options.output_every.unwrap_or(1.0);
        let mut target = output_every;
    
        for i in 1..history.len() - 1 {
            if output_every > 0.0 {
                // print blocks
                if (target - history[i].block).abs() <= (target - history[i + 1].block).abs() {
                    grow::display_history_line(&history[i], &options.output_vars, &options.component, &cycle_indexes);
                    target += output_every;
                } 
            } else if output_every < 0.0 {
                // print cycles
                if i % -output_every as usize == 0 {
                    grow::display_history_line(&history[i], &options.output_vars, &options.component, &cycle_indexes);
                }
            }
        }
    }

    // Display the final moment after failure
    grow::display_history_line(history.last().unwrap(), &options.output_vars, &options.component, &cycle_indexes);
    println!("{}", messages);

    history.to_vec()
}

fn is_input_data_consistent(options: &options::EasiOptions) -> Result<(), String> {
    // Check that k_ut >= material.k1c
    if options.dadn.starts_with("lin_lower_th_5a") ||
        options.dadn.starts_with("lin_lower_th_6a") {
        let k_ut = match options.params.get(&dadn::ParameterLabel::k_ut) {
            Some(value) => value,
            None => {
                error!("Please specify a k_ut value for this dadn");
                std::process::exit(1)
            }
        };
        let k1c = &options.component.material.k1c;

        if k_ut < k1c {
            return Err(format!("k_ut ({}) must not be less than material.k1c ({})", k_ut, k1c));
        }
    }

    if options.rmin >= options.rmax {
        return Err(format!("rmax ({}) must be greater than rmin ({})", options.rmax, options.rmin));
    }

    Ok(())
}

/// Any request for file or info output will result in program
/// termination. This policy is to reduce the complexity for the
/// user as to what the program does.
fn print_and_exit(options: &mut options::EasiOptions) {
    if options.output == TerminatingOutput::List {
        // write out extended list of options and methods
        list::print_list();
        std::process::exit(0);
    }

    if options.scale == 0.0 {
        options.scale = 1.0;
    }
    
    // Process the cycles per the growth method chosen
    let (fatigue_test, unclosed) = options::builder::grow::get_fatigue_test(options);
    let cycles = fatigue_test.get_cycles();

    // Process all the modifications to the sequence
    let mut sequence = cycle::process_seq_mods(&options.sequence, &options.seq_mods);

    // Recreate the sequence after cycle modifications
    if options.seq_mods.cycles {
        sequence.clear();

        for cycle in cycles {
            sequence.push(cycle.min);
            sequence.push(cycle.max);
        }
    }

    // Write out the sequence file.
    if let Some(outfile) = &options.seq_mods.outfile {
        io::write_sequence(outfile, &sequence);
        std::process::exit(0);
    }

    // Write out the cycles file.
    if let Some(outfile) = &options.cycle_mods.outfile {
        io::write_cycles(outfile, cycles);
        std::process::exit(0);
    }

    // write out the beta by converting to a beta table. This can be
    // then read back in using the file: option for beta selection.
    if !options.beta_outfile.is_empty() {
        let geometry = CrackGeometry {
            a: Some(CrackFront{
                length: options.a,
                angle: options.a_angle,
            }),
            c: Some(CrackFront {
                length: options.c,
                angle: options.c_angle,
            }),
            ratio: Some(options.a / options.c),
        };
        let mut beta = beta::get_beta_fn(&options.beta, &options.component, geometry, options.beta_doi.clone());
        let table_beta = beta.as_table();

        // need to write to file
        let path = Path::new(&options.beta_outfile);
        let display = path.display();

        let mut file = match File::create(path) {
            Err(why) => {
                error!(
                    "Error: could not create the file '{}': {}.",
                    display,
                    why
                );
                std::process::exit(1)
            }
            Ok(file) => file,
        };
        let _ = write!(file, "{}", table_beta);
        std::process::exit(0);
    }

    // write out summary information of the sequence
    if options.output == TerminatingOutput::Summary {
        let seq_source = if !options.seq_infile.is_empty() {
            options.seq_infile.clone()
        } else {
            // This is a little vague as the sequence could be either
            // the default sequence or overwritten with a supplied sequence.
            String::from("(Used specified sequence)")
        };
        cycle::summarise_sequence(&seq_source, &sequence, &options.seq_mods);

        let cycle_source = if !options.cycle_infile.is_empty() {
            options.cycle_infile.clone()
        } else {
            format!(
                "(Obtained from sequence using '{:?}' method)",
                options.cycle_method
            )
        };
        cycle::summarise_cycles(
            &cycle_source,
            cycles,
            &unclosed,
            &options.cycle_mods,
        );

        std::process::exit(0)
    }
}
