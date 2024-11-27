use std::collections::{BTreeMap, HashMap};

use fatigue::dadn::ParameterLabel;
use fatigue::grow::OptimisableParameter;
use fatigue::{beta, cycle, fracto, tag};
use crate::options::{EasiOptions, OptimMethod, TerminatingOutput, Verbosity, GrowthMethod};
use crate::optimise::nelder;
use clap::{App, AppSettings, Arg};
use log::error;

use super::VariableAlphaMode;

/// Get the options from the command line.
pub fn get_options_clap(line: &str, options: &mut EasiOptions) {
    let process = App::new("Easigrow: A crack growth modelling program")
        .author("Paul White, Chris Hodgen, Ben Dixon")
        .version(crate_version!())
        .about(include_str!("../description.md"))
        .setting(AppSettings::AllowLeadingHyphen)
        
        .arg(Arg::with_name("method")
             .short("m")
             .long("method")
             .value_name("NAME")
             .possible_values(&GrowthMethod::variants())
             .help("select fatigue growth calculation method (default easigrow)")
             .case_insensitive(true)
             .takes_value(true))

        .arg(Arg::with_name("scale")
             .short("s")
             .long("scale")
             .value_name("SCALE")
             .help("set scale factor to convert load sequence into engineering units. There is no default, so a scale value must be specified by the user.")
             .takes_value(true))
        
        .arg(Arg::with_name("seq_infile")
             .short("q")
             .long("seq_infile")
             .value_name("FILE")
             .help("read a load sequence from file")
             .takes_value(true))
        
        .arg(Arg::with_name("cycle_infile")
             .long("cycle_infile")
             .value_name("FILE")
             .help("read counted cycles from an input file")
             .takes_value(true))

        .arg(Arg::with_name("sequence")
             .long("sequence")
             .value_name("S1,S2,...")
             .help("explicitly define the load sequence at the command line (default: [0.0, 1.0, 0.0])")
             .takes_value(true)
             .require_delimiter(true)
             .allow_hyphen_values(true))

        .arg(Arg::with_name("beta")
             .short("b")
             .long("beta")
             .value_name("NAME")
             .help("select the beta geometry factor model that will be used (default seft-newman84)")
             .takes_value(true))

        .arg(Arg::with_name("beta_doi")
             .long("beta_doi")
             .help("beta's direction of interest. Identify whether a or c is the direction of interest. Only works for beta solutions that allow FCG to be calculated in both the a and c-directions.")
             .possible_values(&["a", "c"])
             .value_name("DIRECTION")
             .takes_value(true))

        .arg(Arg::with_name("beta_outfile")
             .long("beta_outfile")
             .value_name("FILE")
             .help("write a table describing the current beta model to a file")
             .takes_value(true))

        .arg(Arg::with_name("dadn")
             .short("d")
             .long("dadn")
             .value_name("NAME")
             .help("select the dadn vs delta k relationship to be used for a prediction (default white:barter14_aa7050-t7451)")
             .takes_value(true))
                         
        .arg(Arg::with_name("params")
             .short("p")
             .long("parameters")
             .value_name("p1=v1,p2=v2,...,pM=vM")
             .help("set parameters for the chosen da/dN vs. delta k equation")
             .takes_value(true))

        .arg(Arg::with_name("cycle_method")
             .long("cycle_method")
             .help("select method used to count cycles from the load sequence (default rainflow)")
             .possible_values(&["rainflow", "tension"])
             .value_name("METHOD")
             .takes_value(true))

        .arg(Arg::with_name("astart")
             .long("astart")
             .value_name("LENGTH")
             .help("set the initial crack size in the forward (thickness) direction (default 10e-6 m)")
             .takes_value(true))

        .arg(Arg::with_name("cstart")
             .long("cstart")
             .value_name("LENGTH")
             .help("set the initial crack size in the sideways (width) direction (default 10e-6 m)")
             .takes_value(true))

        .arg(Arg::with_name("a_angle")
             .long("a_angle")
             .value_name("ANGLE")
             .help("change the angle associated with crack length in the a-direction (default 0 deg)")
             .takes_value(true))

        .arg(Arg::with_name("c_angle")
             .long("c_angle")
             .value_name("ANGLE")
             .help("change the angle associated with crack length in the c-direction (default 90 deg)")
             .takes_value(true))
        
        .arg(Arg::with_name("forward")
             .long("forward")
             .value_name("DISTANCE")
             .help("set distance from the crack origin to the nearest free edge in the a-direction (thickness) of growth (default: infinite)")
             .takes_value(true))
        
        .arg(Arg::with_name("sideways")
             .long("sideways")
             .value_name("DISTANCE")
             .help("set the distance from the crack origin to nearest free edge in the c-direction (width) of growth (default: infinite)")
             .takes_value(true))
        
        .arg(Arg::with_name("radius")
             .long("radius")
             .value_name("R")
             .help("specify the radius of hole or notch (default: infinite)")
             .takes_value(true))

        .arg(Arg::with_name("rmax")
             .long("rmax")
             .value_name("MAX")
             .help("set maximum allowable R value. Any R greater than this threshold will be set to this value (default: 0.99) Currently only available for: forman, lin_lower_th_5, lin_lower_th_5a, lin_lower_th_5b, lin_lower_th_6, lin_lower_th_6a and file (i.e. tabular data)")
             .takes_value(true))

        .arg(Arg::with_name("rmin")
             .long("rmin")
             .value_name("MIN")
             .help("set minimum allowable R value. Any R less than this threshold will be set to this value (default: -1.0) Currently only available for: forman, lin_lower_th_5, lin_lower_th_5a, lin_lower_th_5b, lin_lower_th_6, lin_lower_th_6a and file (i.e. tabular data)")
             .takes_value(true))

        .arg(Arg::with_name("kneg")
             .long("kneg")
             .help("for R < 0 make delta k = Kmax - Kmin (default: delta k = Kmax for R < 0) Currently only available for: forman, lin_lower_th_5, lin_lower_th_5a, lin_lower_th_5b, lin_lower_th_6, lin_lower_th_6a and file (i.e. tabular data)"))

        .arg(Arg::with_name("kmax_th_zero")
             .long("kmax_th_zero")
             .help("set dadn's kmax threshold to 0 (default deltak_th)"))

        .arg(Arg::with_name("alpha")
             .long("alpha")
             .value_name("ALPHA")
             .help("factor describing the plastic zone state of stress (default: 1.0)")
             .takes_value(true))

        .arg(Arg::with_name("rso")
             .long("rso")
             .value_name("RSO")
             .help("set shut-off ratio used with Willenborg method (default: 3.0)")
             .takes_value(true))
        
        .arg(Arg::with_name("deltak_th")
             .long("deltak_th")
             .value_name("VALUE")
             .help("change the threshold delta k that must be overcome for crack growth to occur (default: 0 MPa sqrt(m))")
             .takes_value(true))
        
        .arg(Arg::with_name("no_yield_limit")
             .long("no_yield_limit")
             .help("do not check failure against the material's yield strength when growing a crack"))

        .arg(Arg::with_name("no_k1c_limit")
             .long("no_k1c_limit")
             .help("do not check failure against the material's fracture toughness when growing a crack"))
        
    // Termination criteria 
    // only required if you want to terminate the crack calculation when this is exceeded
        .arg(Arg::with_name("aend")
             .long("aend")
             .help("set the final crack size in the forward (thickness) direction (default 1e-3 m)")
             .value_name("LENGTH")
             .takes_value(true))

        .arg(Arg::with_name("cend")
             .long("cend")
             .help("set the final crack size in the sideways (width) direction (default 1e-3 m)")
             .value_name("LENGTH")
             .takes_value(true))

        .arg(Arg::with_name("limit_kc")
             .long("limit_kc")
             .value_name("KC")
             .help("set material fracture toughness used for failure checks (default 33 MPa sqrt(m))")
             .takes_value(true))
             
        .arg(Arg::with_name("limit_yield")
             .long("limit_yield")
             .value_name("STRESS")
             .help("set material yield stress (e.g. for failure checks; default: 450 MPa)")
             .takes_value(true))
        
        .arg(Arg::with_name("limit_block")
             .short("N")
             .long("limit_block")
             .value_name("N")
             .help("set maximum number of spectrum blocks before the fatigue calculation terminates (default +500)")
             .takes_value(true))

        .arg(Arg::with_name("youngs_modulus")
             .long("youngs")
             .value_name("MODULUS")
             .help("set Young's modulus (default 71000 MPa)")
             .takes_value(true))

    // Cycle counting
        .arg(Arg::with_name("seq_max")
             .long("seq_max")
             .value_name("MAX")
             .help("limit load sequence points to this maximum value")
             .takes_value(true))

        .arg(Arg::with_name("seq_min")
             .long("seq_min")
             .value_name("MIN")
             .help("limit load sequence points to this minimum value")
             .takes_value(true))

        .arg(Arg::with_name("seq_rem_big")
             .long("seq_rem_big")
             .value_name("BIG")
             .help("remove any load sequence points greater than this threshold")
             .takes_value(true))

        .arg(Arg::with_name("seq_rem_small")
             .long("seq_rem_small")
             .value_name("SMALL")
             .help("remove any load sequence points below this threshold")
             .takes_value(true))

        .arg(Arg::with_name("seq_cycles")
             .long("seq_cyclemods")
             .help("put the fatigue cycles remaining after any cycle modifications into a sequence that can be output via --seq_outfile.")
             .allow_hyphen_values(true))

        .arg(Arg::with_name("seq_outfile")
             .long("seq_outfile")
             .value_name("FILE")
             .help("write load sequence to a file after any modifications")
             .takes_value(true))
                    
        .arg(Arg::with_name("cycle_max")
             .long("cycle_max")
             .value_name("MAX")
             .help("limit the peaks in all cycles to this maximum value")
             .takes_value(true))

        .arg(Arg::with_name("cycle_min")
             .long("cycle_min")
             .value_name("MIN")
             .help("limit the valleys in all cycles to this minimum value")
             .takes_value(true))

        .arg(Arg::with_name("cycle_rem_big")
             .long("cycle_rem_big")
             .value_name("DELTA")
             .help("remove all cycles with a range bigger than DELTA")
             .takes_value(true))

        .arg(Arg::with_name("cycle_rem_small")
             .long("cycle_rem_small")
             .value_name("DELTA")             
             .help("remove all cycles with a range smaller than DELTA")
             .takes_value(true))

        .arg(Arg::with_name("cycle_deadband")
             .long("cycle_deadband")
             .help("remove all fatigue cycles where the peak and valley both fall within the bounds given by <MIN,MAX>")
             .value_name("MIN,MAX")
             .takes_value(true)
             .require_delimiter(true)
             .allow_hyphen_values(true))

        .arg(Arg::with_name("cycle_outfile")
             .long("cycle_outfile")
             .value_name("FILE")
             .help("write counted cycles after any modifications to a file")
             .takes_value(true))
                    
        .arg(Arg::with_name("seq_reorder")
             .short("r") 
             .long("seq_reorder")
             .help("rotate the load sequence so the maximum occurs at its start and end (default: no reordering)"))

        .arg(Arg::with_name("seq_tp")
             .long("seq_tp")
             .help("remove non-turning points from the sequence. By default, this is done prior to cycle counting"))

        .arg(Arg::with_name("variable_alpha")
             .long("variable_alpha")
             .value_name("p1=v1,p2=v2,...,pM=vM")
             .help("choose to use variable alpha and specify the variable alpha parameters: alpha_min, alpha_max, k, x0")
             .takes_value(true))

        .arg(Arg::with_name("variable_alpha_mode")
             .long("variable_alpha_mode")
             .help("select the independent variable used to calculate variable alpha (default LogGrowthRate)")
             .possible_values(&VariableAlphaMode::variants())
             .case_insensitive(true)
             .value_name("MODE")
             .takes_value(true))

        .arg(Arg::with_name("list")
             .short("l")
             .long("list")
             .help("list all available output parameters, cycle counting methods, beta models, da/dN vs. delta k options and calculation methods"))

        .arg(Arg::with_name("output_every")
             .short("n")
             .long("output_every")
             .value_name("N")
             .help("output predicted crack growth history data every N blocks or -N cycles (default 1)")
             .takes_value(true))

        .arg(Arg::with_name("output_cycles")
             .long("output_cycles")
             .value_name("c1,c2,...")
             .help("request predicted crack growth history output at specific cycles in each spectrum block")
             .takes_value(true)
             .require_delimiter(true))

        .arg(Arg::with_name("output_vars")
             .short("o")
             .long("output_vars")
             .value_name("v1,v2,...")
             .help("set variables to be output via a comma-separated list (default block,a,c)")
             .takes_value(true)
             .require_delimiter(true)
             .allow_hyphen_values(true))
                    
        // .arg(Arg::with_name("verbose")
        //      .short("v")
        //      .long("verbose")
        //      .help("print more information"))

        .arg(Arg::with_name("summary")
             .long("summary")
             .help("print a summary of the load sequence and counted cycles, after any modifications.  The sequence modifications and cycle modifications are also summarised"))

        .arg(Arg::with_name("image_outfile")
             .long("image_outfile")
             .value_name("FILE")
             .help("predict the appearance of a crack surface and write to a file")
             .takes_value(true))
                    
        .arg(Arg::with_name("image_bar")
             .long("image_bar")
             .value_name("LENGTH")
             .help("set size of the scale bar displayed for a predicted crack surface image (default 50e-6 m)")
            .takes_value(true))
        
        .arg(Arg::with_name("image_size")
             .long("image_size")
             .value_name("V,H")
             .help("set size of the predicted crack surface image in pixels (default 300x8000)")
             .takes_value(true)
             .require_delimiter(true)
             .allow_hyphen_values(true))

        .arg(Arg::with_name("image_type")
             .long("image_type")
             .possible_values(&["Sem", "Optical"])
             .value_name("TYPE")
             .help("change the type of crack surface image you want to predict (default sem)")
            .takes_value(true))       

        .arg(Arg::with_name("opt_infile") 
             .long("opt_infile")
             .value_name("FILE")
             .help("select optimisation file that references benchmark measurement sets and prediction commands for each set")
             .takes_value(true))
        
        .arg(Arg::with_name("opt_max") 
             .long("opt_max")
             .value_name("N")
             .help("set maximum number of iterations for optimisation (default 100)")
             .takes_value(true))

        .arg(Arg::with_name("opt_sweep")
             .long("opt_sweep")
             .value_name("f1,...,fM")
             .help("apply M scaling factors to each da/dN vs. delta k parameter and method-related parameter during a sweep optimisation")
             .takes_value(true)
             .require_delimiter(true)
             .allow_hyphen_values(true))
   
        .arg(Arg::with_name("opt_tol") 
             .long("opt_tol")
             .value_name("TOL")
             .help("terminate a Nelder Mead optimisation when the change in error between iterations falls below this threshold (default 1e-3)")
             .takes_value(true))

        .arg(Arg::with_name("opt_nelder") 
             .long("opt_nelder")
             .value_name("p1,...,p5")
             .help("change parameters controlling Nelder-Mead optimisation (default step: 0.1, alpha: 1.0, gamma: 2.0, rho: 0.5, sigma: 0.5)")
             .takes_value(true)
             .require_delimiter(true)
             .allow_hyphen_values(true))

        .arg(Arg::with_name("opt_method")
             .long("opt_method")
             .value_name("METHOD")
             .possible_values(&OptimMethod::variants())
             .help("select optimisation method (default Particle)")
             .case_insensitive(true)
             .takes_value(true))
        
        .arg(Arg::with_name("opt_range_dadn")
             .long("opt_range_dadn")
             .value_name("p1=[v1,v2],p2=[v3],...")
             .help("set minimum and maximum allowable values for da/dN vs. delta k parameters during PSO (defaults are per provided in --dadn or --parameters)")
             .takes_value(true))

        .arg(Arg::with_name("opt_range_method")
             .long("opt_range_method")
             .value_name("p1=[v1,v2],p2=[v3],...")
             .help("set minimum and maximum allowable values for method parameters during PSO (default values provided in method-associated inputs)")
             .takes_value(true))

        .arg(Arg::with_name("unbounded")
             .long("unbounded")
             .help("allow parameters to go outside the upper and lower parameter bounds during particle swarm optimisation"))

        .arg(Arg::with_name("swarm_size") 
             .long("swarm_size")
             .value_name("N")
             .help("change number of particles used for PSO (default 40)")
             .takes_value(true))

        .arg(Arg::with_name("opt_linear_crack")
             .long("opt_linear_crack")
             .help("use raw crack lengths for fatigue crack growth rate calculations during optimisation (default: use log crack lengths)"))

        .arg(Arg::with_name("opt_pred_step")
             .long("opt_pred_step")
             .value_name("N")
             .help("keep every Nth cycle in the predicted crack growth histories during an optimisation (default 10)")
             .takes_value(true))

        .arg(Arg::with_name("opt_min_fitness") 
             .long("opt_min_fitness")
             .value_name("MIN")
             .help("change upper threshold of acceptable fitness/error for Particle Swarm optimisation (default 1e-3)")
             .takes_value(true))
            
        .arg(Arg::with_name("crack_infile") 
             .short("c")
             .long("crack_infile")
             .value_name("FILE")
             .help("select file containing a set of crack growth measurements that will be used as a benchmark during optimisation")
             .takes_value(true))

        .arg(Arg::with_name("crack_weight") 
             .long("crack_weight")
             .value_name("WEIGHT")
             .help("change weighting factor applied to the errors calculated when predicting a set of measurements (during an optimisation) (default 1.0)")
             .takes_value(true));

    // turn the commands into matches
    let matches = if line.is_empty() {
        // this will get them directly from the command line
        process.get_matches()
    } else {
        // get them from the string
        // add command word since it skips first argument
        process.get_matches_from(("easigrow ".to_string() + line.trim()).split(' '))
    };

    // basic crack growth options
    if let Ok(method) = value_t!(matches.value_of("method"), GrowthMethod) {
        options.method = method;
    }
    if let Ok(scale) = value_t!(matches, "scale", f64) {
        options.scale = scale;
    }
    if let Some(file) = matches.value_of("seq_infile") {
        options.seq_infile = file.to_string();
    }
    if let Some(file) = matches.value_of("cycle_infile") {
        options.cycle_infile = file.to_string();
    }
    if let Ok(sequence) = values_t!(matches, "sequence", f64) {
        options.sequence = tag::Tag::from(&sequence);
    }
    if let Some(beta) = matches.value_of("beta") {
        options.beta = beta.to_string();
    }
    if let Some(beta_doi) = matches.value_of("beta_doi") {
        options.beta_doi = match beta_doi {
            "a" => Some(beta::DirectionOfInterest::A),
            "c" => Some(beta::DirectionOfInterest::C),
            _ => { error!("Error: invalid beta direction of interest");
                   std::process::exit(2);
            },
        };
    }
    if let Some(beta_outfile) = matches.value_of("beta_outfile") {
        options.beta_outfile = beta_outfile.to_string();
        options.output = TerminatingOutput::BetaTable;
    }
    if let Some(dadn) = matches.value_of("dadn") {
        options.dadn = dadn.to_string();
    }
    if let Some(params) = matches.value_of("params") {
        options.params = parse_parameters(params);
    }

    if let Some(cycle_method) = matches.value_of("cycle_method") {
        options.cycle_method = match cycle_method {
            "rainflow" => cycle::CycleMethod::Rainflow,
            "tension" => cycle::CycleMethod::Tension,
            _ => { error!("Error: unknown cycle counting method");
                   std::process::exit(2);
            },

        };
    }
    if let Ok(alpha) = value_t!(matches, "alpha", f64) {
        options.alpha = alpha;
    }
    if let Ok(r_so) = value_t!(matches, "rso", f64) {
        options.r_so = r_so;
    }
    if let Ok(deltak_th) = value_t!(matches, "deltak_th", f64) {
        options.deltak_th = Some(deltak_th);
    }
    if let Some(variable_alpha) = matches.value_of("variable_alpha") {
        options.variable_alpha = Some(parse_list_of_str_f64_pairs(variable_alpha));
    }
    if let Ok(variable_alpha_mode) = value_t!(matches.value_of("variable_alpha_mode"), VariableAlphaMode) {
        options.variable_alpha_mode = variable_alpha_mode;
    }

    // crack and geometry options
    if let Ok(astart) = value_t!(matches, "astart", f64) {
        options.a = astart;
    }
    if let Ok(cstart) = value_t!(matches, "cstart", f64) {
        options.c = cstart;
    }
    if let Ok(a_angle) = value_t!(matches, "a_angle", f64) {
        options.a_angle = a_angle * (std::f64::consts::PI / 180.0);
    }
    if let Ok(c_angle) = value_t!(matches, "c_angle", f64) {
        options.c_angle = c_angle * (std::f64::consts::PI / 180.0);
    }
    if let Ok(limit_a) = value_t!(matches, "aend", f64) {
        options.a_limit = limit_a;
    }
    if let Ok(limit_c) = value_t!(matches, "cend", f64) {
        options.c_limit = limit_c;
    }
    if let Ok(forward) = value_t!(matches, "forward", f64) {
        options.component.forward = forward;
    }
    if let Ok(sideways) = value_t!(matches, "sideways", f64) {
        options.component.sideways = sideways;
    }
    if let Ok(radius) = value_t!(matches, "radius", f64) {
        options.component.radius = radius;
    }
    if let Ok(limit_kc) = value_t!(matches, "limit_kc", f64) {
        options.component.material.k1c = limit_kc;
    }
    if let Ok(limit_yield) = value_t!(matches, "limit_yield", f64) {
        options.component.material.yield_stress = limit_yield;
    }
    if let Ok(limit_block) = value_t!(matches, "limit_block", f64) {
        options.block_limit = limit_block;
    }
    if let Ok(youngs_modulus) = value_t!(matches, "youngs_modulus", f64) {
        options.component.material.youngs_modulus = youngs_modulus;
    }
    if let Ok(rmax) = value_t!(matches, "rmax", f64) {
        options.rmax = rmax;
    }
    if let Ok(rmin) = value_t!(matches, "rmin", f64) {
        options.rmin = rmin;
    }
    if matches.is_present("kneg") {
        options.kneg = true;
    }
    if matches.is_present("no_yield_limit") {
        options.termination_tests.limit_yield = false;
    }
    if matches.is_present("no_k1c_limit") {
        options.termination_tests.limit_k1c = false;
    }
    if matches.is_present("kmax_th_zero") {
        options.kmax_th_zero = true;
    }

    // ssequence modifications
    if let Ok(max) = value_t!(matches, "seq_max", f64) {
        options.seq_mods.cap_max = Some(max);
    }
    if let Ok(min) = value_t!(matches, "seq_min", f64) {
        options.seq_mods.cap_min = Some(min);
    }
    if let Ok(big) = value_t!(matches, "seq_rem_big", f64) {
        options.seq_mods.remove_bigger = Some(big);
    }
    if let Ok(small) = value_t!(matches, "seq_rem_small", f64) {
        options.seq_mods.remove_smaller = Some(small);
    }
    if matches.is_present("seq_cycles") {
        options.seq_mods.cycles = true;
    }
    if let Some(output) = matches.value_of("seq_outfile") {
        options.seq_mods.outfile = Some(output.to_string());
        options.output = TerminatingOutput::SequenceFile;
    }
    if matches.is_present("seq_reorder") {
        options.seq_mods.reorder = true;
    }
    if matches.is_present("seq_tp") {
        options.seq_mods.turning_points = true;
    }

    // cycle modifications
    if let Ok(max) = value_t!(matches, "cycle_max", f64) {
        options.cycle_mods.cap_max = Some(max);
    }
    if let Ok(min) = value_t!(matches, "cycle_min", f64) {
        options.cycle_mods.cap_min = Some(min);
    }
    if let Ok(big) = value_t!(matches, "cycle_rem_big", f64) {
        options.cycle_mods.remove_bigger = Some(big);
    }
    if let Ok(small) = value_t!(matches, "cycle_rem_small", f64) {
        options.cycle_mods.remove_smaller = Some(small);
    }
    if let Ok(deadband) = values_t!(matches, "cycle_deadband", f64) {
        options.cycle_mods.remove_region = Some(cycle::Limit::from_vec(&deadband));
    }
    if let Some(output) = matches.value_of("cycle_outfile") {
        options.cycle_mods.outfile = Some(output.to_string());
        options.output = TerminatingOutput::CycleFile;
    }

    // output options
    if let Ok(output_every) = value_t!(matches, "output_every", f64) {
        options.output_every = Some(output_every);
    }
    if let Ok(output_vars) = values_t!(matches, "output_vars", String) {
        options.output_vars = output_vars;
    }
    if let Ok(output_cycles) = values_t!(matches, "output_cycles", usize) {
        options.output_cycles = output_cycles;
    }
    if matches.is_present("list") {
        options.output = TerminatingOutput::List;
    }
    if matches.is_present("summary") {
        options.output = TerminatingOutput::Summary;
    }
    if matches.is_present("verbose") {
        options.verbosity = Verbosity::Verbose;
        println!("Matches: {:?}", matches);
    }

    // Image generation options.
    if let Some(file) = matches.value_of("image_outfile") {
        options.image.file = file.to_string();
    }
    if let Ok(image_bar) = value_t!(matches, "image_bar", f64) {
        options.image.barlength = image_bar;
    }
    if matches.is_present("image_type") {
        options.image.image =
            value_t!(matches, "image_type", fracto::ImageType).unwrap_or_else(|e| e.exit())
    }
    if let Ok(image_size) = values_t!(matches, "image_size", u32) {
        options.image.xsize = image_size[0];
        options.image.ysize = image_size[1];
    }

    // Optimisation Options
    if let Some(file) = matches.value_of("opt_infile") {
        options.optimise.file = file.to_string();
    }
    if let Ok(opt_method) = value_t!(matches.value_of("opt_method"), OptimMethod) {
        options.optimise.method = opt_method;
    };
    if let Ok(opt_sweep) = values_t!(matches, "opt_sweep", f64) {
        options.optimise.sweep = opt_sweep;
    }
    if let Ok(opt_nelder) = values_t!(matches, "opt_nelder", f64) {
        options.optimise.nelder_params = nelder::Parameters::new(&opt_nelder);
    }
    if let Ok(opt_max) = value_t!(matches, "opt_max", usize) {
        options.optimise.maxiter = opt_max;
    }
    if let Ok(opt_tol) = value_t!(matches, "opt_tol", f64) {
        options.optimise.tol = opt_tol;
    }
    if let Some(opt_range_dadn) = matches.value_of("opt_range_dadn") {
        options.optimise.range_dadn = parse_range_dadn(opt_range_dadn);
    }
    if let Some(opt_range_method) = matches.value_of("opt_range_method") {
        options.optimise.range_method = parse_range_method(opt_range_method);
    }
    if matches.is_present("unbounded") {
        options.optimise.bounds_constraints = false;
    }
    if let Ok(swarm_size) = value_t!(matches, "swarm_size", usize) {
        options.optimise.swarm_size = swarm_size;
    }
    if matches.is_present("opt_linear_crack") {
        options.optimise.use_log_crack_depth = false;
    }
    if let Ok(opt_keep_every) = value_t!(matches, "opt_pred_step", u32) {
        options.optimise.keep_every = opt_keep_every;
    }
    if let Ok(opt_min_fitness) = value_t!(matches, "opt_min_fitness", f64) {
        options.optimise.min_fitness = opt_min_fitness;
    }

    // Fracto options
    if let Some(file) = matches.value_of("crack_infile") {
        options.crack_infile = file.to_string();
    }
    if let Ok(crack_weight) = value_t!(matches, "crack_weight", f64) {
        options.crack_weight = crack_weight;
    }
}


// Input expected as a comma separated list of paired values, of the 
// form str=f64
// e.g a=1.0,b=1e-6,c=-0.4
fn parse_parameters(input: &str) -> BTreeMap<ParameterLabel, f64> {
    let mut result = BTreeMap::new();
    // Split by comma to get vector of paired values
    let pairs = input.split(',').collect::<Vec<&str>>();

    // Iterate over each pair and insert into map
    for pair in pairs {
        let key_value = pair.split('=').collect::<Vec<&str>>();
        if key_value.len() != 2 {
            error!("Error: Invalid parameter format: {}", pair);
            std::process::exit(1)
        }

        let key = match ParameterLabel::from_text(key_value[0]) {
            Some(label) => {
                label
            },
            None => {
                error!("Error: Unknown parameter label: {}", key_value[0]);
                std::process::exit(1)
            }
        };

        let value = match key_value[1].parse() {
            Ok(v) => v,
            Err(why) => {
                error!("Error: Invalid parameter value: {}. {}", key_value[1], why);
                std::process::exit(1)
            }
        }; 

        result.insert(key, value);
    }

    result
}


// Input expected as a comma separated list of paired values, of the 
// form str=f64
// e.g a=1.0,b=1e-6,c=-0.4
fn parse_list_of_str_f64_pairs(input: &str) -> HashMap<String, f64> {
    let mut result = HashMap::new();

    // Split by comma to get vector of paired values
    let pairs = input.split(',').collect::<Vec<&str>>();

    // Iterate over each pair and insert into map
    for pair in pairs {
        let key_value = pair.split('=').collect::<Vec<&str>>();
        if key_value.len() != 2 {
            error!("Error: Invalid parameter format: {}", pair);
            std::process::exit(1)
        }

        let key = key_value[0].to_string();

        let value = match key_value[1].parse() {
            Ok(v) => v,
            Err(why) => {
                error!("Error: Invalid parameter value: {}. {}", key_value[1], why);
                std::process::exit(1)
            }
        };

        result.insert(key, value);
    }

    result
}


fn parse_range_dadn(input: &str) -> HashMap<ParameterLabel, (f64, f64)> {
    let mut result = HashMap::new();

    let first_pass = parse_opt_range(input);

    for (k, v) in first_pass {
        let key = match ParameterLabel::from_text(&k) {
            Some(label) => {
                label
            },
            None => {
                error!("Error: Unknown parameter label: {}", k);
                std::process::exit(1)
            }
        };

        result.insert(key, v);
    }

    result
}


fn parse_range_method(input: &str) -> HashMap<OptimisableParameter, (f64, f64)> {
    let mut result = HashMap::new();

    let first_pass = parse_opt_range(input);

    for (k, v) in first_pass {
        let key = match OptimisableParameter::from_text(&k) {
            Some(label) => {
                label
            },
            None => {
                error!("Error: Unknown parameter label: {}", k);
                std::process::exit(1)
            }
        };

        result.insert(key, v);
    }

    result
}

// Expected format:
// a=[1,2],b=[1.0,2.0],c=[-1]
fn parse_opt_range(input: &str) -> HashMap<String, (f64, f64)> {
    let mut result = HashMap::new();

    let parse = |input: &str| -> f64 {
        match input.parse() {
            Ok(value) => value,
            Err(why) => {
                error!("Error: Invalid range parameter value: {}. {}", input, why);
                std::process::exit(1)
            }
        }
    };

    let mut building_key = true;
    let mut building_values = false;
    let mut key = "".to_string();
    let mut value = "".to_string();
    let mut values = vec![];
    for ch in input.chars() {
        match ch {
            '=' => building_key = false,
            '[' => building_values = true,
            ']' => {
                values.push(parse(&value));
                value.clear();

                let result_values = match values.len() {
                    1 => (values[0], values[0]),
                    2 => (values[0], values[1]),
                    _ => {
                        error!("Error: Incorrect number of range values. Expected 1 or 2, received {}", values.len());
                        std::process::exit(1);
                    } 
                };
                result.insert(key.to_string(), result_values);

                key.clear();
                values.clear();
                building_values = false;
            }
            ',' => {
                if building_values {
                    values.push(parse(&value));
                    value.clear()
                } else {
                    building_key = true
                }
            }
            _ if building_key => key.push(ch),
            _ if building_values => value.push(ch),
            _ => ()
        }
    }

    result
}


#[cfg(test)]
mod tests {
    use super::parse_opt_range;

    
    #[test]
    fn parse_opt_range_works_with_single_value() {
        let input = "a=[1]";

        let result = parse_opt_range(input);
        let values = result.get("a").unwrap();

        assert!(values.0 - 1.0 <= f64::EPSILON);
        assert!(values.1 - 1.0 <= f64::EPSILON);
    }

    #[test]
    fn parse_opt_range_works_on_single_parameter() {
        let input = "a=[1,2]";

        let result = parse_opt_range(input);
        let values = result.get("a").unwrap();
        
        assert!(values.0 - 1.0 <= f64::EPSILON);
        assert!(values.1 - 2.0 <= f64::EPSILON);
    }

    #[test]
    fn parse_opt_range_works_on_multiple_parameters() {
        let input = "a=[1,2],b=[-3,-4]";

        let result = parse_opt_range(input);
        let values = result.get("a").unwrap();
        
        assert!(values.0 - 1.0 <= f64::EPSILON);
        assert!(values.1 - 2.0 <= f64::EPSILON);

        let values = result.get("b").unwrap();
        
        assert!(values.0 - -3.0 <= f64::EPSILON);
        assert!(values.1 - -4.0 <= f64::EPSILON);
    }
}
