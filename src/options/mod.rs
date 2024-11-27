/// These are the data structures for command line options as well as
/// the initial default values.

use std::collections::{BTreeMap, HashMap};
use std::f64;
use std::f64::consts::FRAC_PI_2;
use std::string::String;
use fatigue::dadn::{self, ParameterLabel};
use fatigue::grow::OptimisableParameter;
use fatigue::{beta, cycle, fracto, grow, io, material, tag};
use fatigue::COMMENT;

use crate::optimise::nelder;

pub mod builder;
pub mod clap;

arg_enum!{
    #[derive(Debug, Clone, PartialEq)]
    pub enum OptimMethod {
        Sweep,
        Nelder,
        Particle,
    }
}

arg_enum!{
    #[derive(Debug, Clone)]
    pub enum CycleMethod {
        Rainflow,
        Tension
    }
}

arg_enum!{
    #[derive(Debug, Clone)]
    pub enum GrowthMethod {
        Easigrow,
        Willenborg
    }
}

arg_enum! {
    #[derive(Debug, Clone)]
    pub enum VariableAlphaMode {
        GrowthRate,
        LogGrowthRate,
        CrackLength,
        LogCrackLength,
    }
}

#[derive(Debug, Clone)]
pub struct Optimise {
    pub file: String,
    pub method: OptimMethod,
    pub maxiter: usize,
    pub tol: f64,
    pub nelder_params: nelder::Parameters,
    pub sweep: Vec<f64>,
    pub range_dadn: HashMap<ParameterLabel, (f64, f64)>,
    pub range_method: HashMap<OptimisableParameter, (f64, f64)>,
    pub swarm_size: usize,
    pub bounds_constraints: bool,
    pub use_log_crack_depth: bool,
    /// Keep every Nth history entry during optimisation
    pub keep_every: u32,
    pub min_fitness: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub enum Verbosity {
    Verbose,
    Terse,
}
#[derive(Debug, Clone, PartialEq)]
pub enum TerminatingOutput {
    List,
    Summary,
    BetaTable,
    SequenceFile,
    CycleFile,
    None,
}

/// Option data for performing a crack growth calculation.
#[derive(Debug, Clone)]
pub struct EasiOptions {
    /// Specify which growth method will be used
    pub method: GrowthMethod,
    /// Initial starting size in the forward direction
    pub a: f64,
    /// Initial starting size in the sideways direction
    pub c: f64,
    /// Angle of the crack in the forward direction
    pub a_angle: f64,
    /// Angle of the crack in the sideways direction
    pub c_angle: f64,
    /// Final crack depth in the forward direction
    pub a_limit: f64,
    /// Final crack depth in the sideways direction
    pub c_limit: f64,
    /// Maximum number of cycles to run for.
    pub block_limit: f64,
    /// Type of output
    pub output: TerminatingOutput,
    /// Level of verbosity
    pub verbosity: Verbosity,
    /// Name of inbuilt dadn data to use.
    pub dadn: String,
    /// Dadn parameters to be used in optimisation or any crack calculation.
    // A BTreeMap is being used to guarantee ordering during iteration for
    // use in optimisations. 
    pub params: BTreeMap<ParameterLabel, f64>,
    /// Parameters to be written out.
    pub output_vars: Vec<String>,
    /// Output 'every' block.
    pub output_every: Option<f64>,
    /// Cycle positions to write out
    pub output_cycles: Vec<usize>,
    /// Name of damage model to use.
    pub cycle_method: cycle::CycleMethod,
    /// Scale the sequence by this factor, typically stress
    pub scale: f64,
    /// Beta model to use.
    pub beta: String,
    pub beta_outfile: String,
    /// Information describing the shape of the component.
    pub component: grow::Component,
    /// Sequence to be used.
    pub sequence: Vec<tag::Tag>,
    /// Simplified version of the measured crack file.
    pub fracto: Vec<io::Measurement>,
    /// Data from crackfile.
    pub cycles: Vec<cycle::Cycle<tag::Tag>>,
    /// Info for optimisation.
    pub optimise: Optimise,
    /// Data for generating a reconstructed image.
    pub image: fracto::ImageData,
    /// Name of crack growth data file used for target in optimisation.
    pub crack_infile: String,
    /// Weighting factor for crack errors
    pub crack_weight: f64,
    /// Sequence information.
    pub seq_mods: cycle::SequenceModifiers,
    pub cycle_mods: cycle::CycleModifiers,
    /// Filename of sequence.
    pub seq_infile: String,
    /// Name of cycle file for inputting data.
    pub cycle_infile: String,
    /// Maximum R (kmin/kmax) value
    pub rmax: f64,
    /// Minimum R (kmin/kmax) value
    pub rmin: f64,
    /// When true, use negative kmin values when calculating delta k
    pub kneg: bool,
    /// State of stress used during plastic zone calculations
    pub alpha: f64,
    /// Shut off ratio used with the Willenborg method
    pub r_so: f64,
    /// Delta k threshold for use with Willenborg method
    /// This option takes lower priority in cases where a dadn equation specifies a deltak_th parameter
    pub deltak_th: Option<f64>,
    /// Control the ability to run different growth termination tests
    pub termination_tests: grow::TerminationTests,
    /// Variable alpha parameters
    pub variable_alpha: Option<HashMap<String, f64>>,
    /// Which variable will be used in variable alpha calculations
    pub variable_alpha_mode: VariableAlphaMode,
    /// Flag to set dadn's kmax threshold to 0, rather than deltak_th
    pub kmax_th_zero: bool,
    /// Direction of interest for some betas can be changed
    pub beta_doi: Option<beta::DirectionOfInterest>,
}

impl EasiOptions {
    /// Returns true if dadn is from a file, false in all other cases
    pub fn is_dadn_from_file(&self) -> bool {
        self.dadn.starts_with("file")
    }
}

impl std::fmt::Display for EasiOptions {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let _e = writeln!(f, "{}method: {:?}", COMMENT, self.method);
        let _e = writeln!(f, "{}a: {:?}", COMMENT, self.a);
        let _e = writeln!(f, "{}a_limit: {:?}", COMMENT, self.a_limit);
        let _e = writeln!(f, "{}block_limit: {:?}", COMMENT, self.block_limit);
        let _e = writeln!(f, "{}params: {:?}", COMMENT, self.params);
        let _e = writeln!(f, "{}output_vars: {:?}", COMMENT, self.output_vars);
        let _e = writeln!(f, "{}output_every: {:?}", COMMENT, self.output_every);
        let _e = writeln!(f, "{}output_cycles: {:?}", COMMENT, self.output_cycles);
        let _e = writeln!(f, "{}scale: {:?}", COMMENT, self.scale);
        let _e = writeln!(f, "{}beta: {:?}", COMMENT, self.beta);
        let _e = writeln!(f, "{}component: {:?}", COMMENT, self.component);
        let _e = writeln!(f, "{}seq_infile: {:?}", COMMENT, self.seq_infile);
        let _e = writeln!(f, "{}seq_mods: {:?}", COMMENT, self.seq_mods);
        let _e = writeln!(f, "{}cycle_method: {:?}", COMMENT, self.cycle_method);
        let _e = writeln!(f, "{}cycle_infile: {:?}", COMMENT, self.cycle_infile);
        let _e = writeln!(f, "{}cycle_mods: {:?}", COMMENT, self.cycle_mods);

        if !self.optimise.file.is_empty() {
            let _e = writeln!(f, "{}optimise: {:?}", COMMENT, self.optimise);
            let _e = writeln!(f, "{}crack_infile: {:?}", COMMENT, self.crack_infile);
            let _e = writeln!(f, "{}crack_weight: {:?}", COMMENT, self.crack_weight);
        }

        if !self.image.file.is_empty() {
            let _e = writeln!(f, "{}fracto: {:?}", COMMENT, self.fracto);
            let _e = writeln!(f, "{}image: {:?}", COMMENT, self.image);
        }
        write!(f, "{}dadn: {:?}", COMMENT, self.dadn)
    }
}

/// read in each sequence and the measured crack growth file
///
/// This stuff only needs to be performed once such as populating the
/// values of tables read from files.
pub fn read_all_files(options: &mut EasiOptions) {
    if !options.seq_infile.is_empty() {
        options.sequence = io::read_sequence(&options.seq_infile);
    }

    if !options.cycle_infile.is_empty() {
        options.cycles = io::read_afgrow_cycles(&options.cycle_infile);
    }

    if !options.crack_infile.is_empty() {
        options.fracto = io::read_fracto_file(&options.crack_infile, options.sequence.len());
    }
}

pub fn get_default_options() -> EasiOptions {
    // default options for easigrow
    EasiOptions {
        // crack growth formulation
        method: GrowthMethod::Easigrow,
        a: 10e-6,
        c: 10e-6,
        a_angle: 0.0,
        c_angle: FRAC_PI_2,
        dadn: "white:barter14-aa7050t7451".to_string(),
        params: BTreeMap::new(),
        beta: "seft-newman84".to_string(),
        beta_doi: None,
        beta_outfile: "".to_string(),
        component: grow::Component {
            sideways: f64::INFINITY,
            forward: f64::INFINITY,
            radius: f64::INFINITY,
            material: material::Properties::default(),
        },
        rmax: dadn::Options::default().rmax,
        rmin: dadn::Options::default().rmin,
        kneg: dadn::Options::default().kneg,
        alpha: 1.0,
        r_so: 3.0,
        deltak_th: None,
        variable_alpha: None,
        variable_alpha_mode: VariableAlphaMode::LogGrowthRate,
        kmax_th_zero: false,

        // termination criteria
        a_limit: 1e-3,
        c_limit: 1e-3, // (m)
        block_limit: 500.0,
        termination_tests: grow::TerminationTests {
            limit_k1c: true,
            limit_yield: true,
        },

        // sequence info
        scale: 0.0,                                 // MPa
        seq_infile: "".to_string(),                 // name of sequence file
        sequence: tag::Tag::from(&[0.0, 1.0, 0.0]), // Constant amplitude cycle
        cycles: vec![],                             // cycles from
        cycle_infile: "".to_string(),
        cycle_method: cycle::CycleMethod::Rainflow,
        seq_mods: cycle::SequenceModifiers {
            cap_max: None,
            cap_min: None,
            remove_bigger: None,
            remove_smaller: None,
            cycles: false,
            reorder: false,
            turning_points: false,
            outfile: None,
        },
        cycle_mods: cycle::CycleModifiers {
            cap_max: None,
            cap_min: None,
            remove_bigger: None,
            remove_smaller: None,
            remove_region: None,
            outfile: None,
        },

        // crack comparison
        crack_infile: "".to_string(), // name of fracto file
        crack_weight: 1.0,
        fracto: vec![],

        // solution type
        optimise: Optimise {
            file: "".to_string(),
            maxiter: 100,
            method: OptimMethod::Particle,
            tol: 1e-4,
            nelder_params: nelder::Parameters::default(),
            sweep: vec![0.8, 1.0, 1.25],
            // bound_min: BTreeMap::new(),
            // bound_max: BTreeMap::new(),
            range_dadn: HashMap::new(),
            range_method: HashMap::new(),
            swarm_size: 40,
            bounds_constraints: true,
            use_log_crack_depth: true,
            keep_every: 10,
            min_fitness: 1e-3,
        },

        image: fracto::ImageData {
            file: "".to_string(),
            barlength: 50e-6,
            xsize: 300,
            ysize: 8000,
            image: fracto::ImageType::Sem,
        },

        // output options
        output_vars: vec!["block".to_string(), "a".to_string(), "c".to_string()],
        output_every: None,
        output_cycles: vec![],
        output: TerminatingOutput::None,
        verbosity: Verbosity::Terse,
    }
}
