//! Grow a fatigue crack from initial crack size until failure.
#![allow(clippy::useless_let_if_seq)]

use crate::beta::{BetaResult, DirectionOfInterest};
use crate::cycle;
use crate::cycle::Cycle;
use crate::dadn::ParameterLabel;
use crate::material;
use crate::tag;
use crate::tag::Tag;
use crate::COMMENT;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::f64::consts::PI;
use std::fmt::{self, Write};
use std::process;

pub mod easigrow;
pub mod stress_constraint;
pub mod willenborg;

/// Common interface for the various crack growth methods
pub trait FatigueTest {
    /// Run a full growth prediction
    fn run(&mut self) -> (&[History], String);

    /// Produce the minimal required data for optimisation
    fn run_for_optimisation(&mut self) -> &[BasicHistory];

    fn update_parameters(&mut self, parameters: &[f64]);

    fn reset(&mut self);

    fn get_component(&self) -> &Component;

    fn get_cycles(&self) -> &[Cycle<Tag>];

    fn get_history(&self) -> &[History];

    fn inner_clone(&self) -> Box<dyn FatigueTest + Send + Sync>;

    fn printable_parameters(&self, parameters: &[f64]) -> String;
}

impl Clone for Box<dyn FatigueTest + Send + Sync> {
    fn clone(&self) -> Self {
        self.inner_clone()
    }
}

/// Data collected for each crack growth prediction cycle.
#[derive(Debug, Clone)]
pub struct History {
    /// Floating point number for block where fractional parts gives fraction of cycle in the block
    pub block: f64,
    /// Applied scaling stress for this cycle.
    pub stress: f64,
    /// Cycle information.
    pub cycle: Cycle<tag::Tag>,
    /// Stress intensity around the crack
    pub kmax: StressIntensity,
    pub kmin: StressIntensity,
    /// stress intensity range
    pub dk: StressIntensityDelta,
    /// beta values around the crack front
    pub beta: BetaResult,
    /// growth increment around the crack front
    pub length_deltas: CrackLengthDelta,
    /// overall geometry
    pub geometry: CrackGeometry,
}

/// A stripped down version of history to be used with optimisations
#[derive(Debug, Clone)]
pub struct BasicHistory {
    pub block: f64,
    pub crack_length: f64,
}

/// Used for storing various types of lengths
#[derive(Debug, Clone)]
pub struct CrackLength {
    pub a: Option<f64>,
    pub c: Option<f64>,
}

// The aliases below exist to make code clearer. It makes more sense
// to do this instead of duplicating the exact same struct layout
// (and implementations) multiple times.

/// Useful for representing deltas specifically
pub type CrackLengthDelta = CrackLength;

/// Represent a stress intensity (kmin, kmax)
pub type StressIntensity = CrackLength;

/// Represent a stress intensity delta (delta k)
pub type StressIntensityDelta = CrackLength;

/// A single crack front
#[derive(Debug, Clone)]
pub struct CrackFront {
    pub length: f64,
    pub angle: f64,
}

/// Due to varying geometries/betas, there are different types
/// of crack front combinations available
#[derive(Debug, Clone)]
pub struct CrackGeometry {
    /// Crack front in the `forward` or `d` direction
    pub a: Option<CrackFront>,
    /// Crack front in the `sideways` or `b` or `w` direction
    pub c: Option<CrackFront>,
    /// Used in some beta functions where an a\c ratio is fixed
    pub ratio: Option<f64>,
}

#[derive(Debug, Clone, Default)]
/// describes the geometry and material containing the crack
pub struct Component {
    pub forward: f64,
    pub sideways: f64,
    pub radius: f64,
    pub material: material::Properties,
}

#[derive(Debug, Clone, Default)]
pub struct TerminationTests {
    pub limit_yield: bool,
    pub limit_k1c: bool,
}

/// Common data required by all methods
#[derive(Debug, Clone, Default)]
pub struct Grow {
    pub cycles: Vec<Cycle<Tag>>,
    pub component: Component,
    pub scale: f64,
    pub limit: CrackLength,
    pub block_limit: f64,
    pub termination_tests: TerminationTests,
    pub keep_every: u32,
}

/// A common collection of parameters used by growth methods which can be optimised
#[allow(non_camel_case_types)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum OptimisableParameter {
    rso,
    alpha,
    deltak_th,
    alpha_min,
    alpha_max,
    k,
    x0,
}


impl OptimisableParameter {
    pub fn text(&self) -> &'static str {
        match *self {
            OptimisableParameter::rso => "rso",
            OptimisableParameter::alpha => "alpha",
            OptimisableParameter::deltak_th => "deltak_th",
            OptimisableParameter::alpha_min => "alpha_min",
            OptimisableParameter::alpha_max => "alpha_max",
            OptimisableParameter::k => "k",
            OptimisableParameter::x0 => "x0",
        }
    }

    pub fn from_text(input: &str) -> Option<OptimisableParameter> {
        match input {
            "rso" => Some(OptimisableParameter::rso),
            "alpha" => Some(OptimisableParameter::alpha),
            "deltak_th" => Some(OptimisableParameter::deltak_th),
            "alpha_min" => Some(OptimisableParameter::alpha_min),
            "alpha_max" => Some(OptimisableParameter::alpha_max),
            "k" => Some(OptimisableParameter::k),
            "x0" => Some(OptimisableParameter::x0),
            _ => None,
        }
    }
}


//TODO: Implementation
impl fmt::Debug for dyn FatigueTest + Send + Sync {
    fn fmt(&self, _f: &mut fmt::Formatter) -> fmt::Result {
        todo!()
    }
}

impl History {
    // Get an initial history entry
    pub fn new(beta: BetaResult, geometry: CrackGeometry) -> History {
        History {
            block: 0.0,
            stress: 0.0,
            cycle: cycle::Cycle {
                max: tag::Tag {
                    value: 0.0,
                    index: usize::MAX,
                },
                min: tag::Tag {
                    value: 0.0,
                    index: usize::MAX,
                },
            },
            kmax: StressIntensity::default(),
            kmin: StressIntensity::default(),
            dk: StressIntensityDelta::default(),
            beta,
            length_deltas: CrackLengthDelta::default(),
            geometry,
        }
    }
}

impl CrackLength {
    pub fn new(a: Option<f64>, c: Option<f64>) -> Self {
        Self { a, c }
    }
}

impl Default for CrackLength {
    fn default() -> Self {
        Self::new(None, None)
    }
}

impl CrackGeometry {
    pub fn get_length_a(&self) -> Option<f64> {
        self.a.as_ref().map(|a| a.length)
    }

    pub fn get_length_c(&self) -> Option<f64> {
        self.c.as_ref().map(|c| c.length)
    }

    /// Update the crack length in the `a` direction
    ///
    /// If there is no `a` direction set, this method will have no effect
    pub fn set_length_a(&mut self, length: f64) {
        if self.a.is_some() {
            self.a.as_mut().unwrap().length = length;
        }
    }

    /// Update the crack length in the `c` direction
    ///
    /// If there is no `c` direction set, this method will have no effect
    pub fn set_length_c(&mut self, length: f64) {
        if self.c.is_some() {
            self.c.as_mut().unwrap().length = length;
        }
    }
}

pub fn k_on_stress(beta: f64, crackfront_length: f64) -> f64 {
    beta * (PI * crackfront_length).sqrt()
}

pub fn k(stress: f64, beta: f64, crackfront_length: f64) -> f64 {
    stress * k_on_stress(beta, crackfront_length)
}

pub fn kmin(cycle: &Cycle<Tag>, beta: f64, crackfront_length: f64) -> f64 {
    cycle.min.value * k_on_stress(beta, crackfront_length)
}

pub fn kmax(cycle: &Cycle<Tag>, beta: f64, crackfront_length: f64) -> f64 {
    cycle.max.value * k_on_stress(beta, crackfront_length)
}

/// Extract the length from a given geometry in the direction of interest
pub fn length_of_interest(
    direction_of_interest: &DirectionOfInterest,
    geometry: &CrackGeometry,
) -> f64 {
    match direction_of_interest {
        DirectionOfInterest::A => geometry.a.as_ref().unwrap().length,
        DirectionOfInterest::C => geometry.c.as_ref().unwrap().length,
    }
}

pub fn display_history_header(output: &[String]) {
    // write out the headers
    if !output.is_empty() {
        for out in output.iter() {
            print!("{:>12} ", out);
        }
        println!();
    }
}

/// Test to see if a history line should be output.
pub fn output_cycle_history(
    his: &History,
    previous_block: f64,
    every: i32,
    output_lines: &BTreeSet<usize>,
    cycle_no: usize,
) -> bool {
    // if every is positive write out every nth block, otherwise if
    // every is negative write out every nth cycle
    // if every == 0, then print out each block with a whole, integer part
    // greater than the previous block. i.e. the first cycle of each block
    let frequency = (every > 0 && his.block as i32 % every == 0)
        || (every < 0 && cycle_no % -every as usize == 0)
        || (every == 0 && his.block.floor() > previous_block.floor());

    // output only if the cycle constains the specific sequence line
    if !output_lines.is_empty() && every > 0 {
        frequency
            && (output_lines.contains(&his.cycle.max.index)
                || output_lines.contains(&his.cycle.min.index))
    } else {
        frequency
    }
}

// print a line of the history data
pub fn display_history_line(his: &History, output: &[String], component: &Component, cycle_indexes: &HashMap<usize, usize>) {
    for out in output {
        let a_on_c = match his.geometry.ratio {
            Some(ratio) => Some(ratio),
            None => {
                if his.geometry.a.is_some() && his.geometry.c.is_some() {
                    Some(
                        his.geometry.a.as_ref().unwrap().length
                            / his.geometry.c.as_ref().unwrap().length,
                    )
                } else {
                    None
                }
            }
        };

        match out.trim() {
            "loadline_no" => if his.cycle.max.index == usize::MAX {
                print!("{:>12} ", "N/A")
            } else {
                print!("{:12} ", his.cycle.max.index)
            },
            "cycle_no" => match cycle_indexes.get(&his.cycle.max.index) {
                Some(cycle_no) => print!("{:12.4} ", cycle_no),
                None => print!("{:>12} ", "N/A"),
            },
            "block" => print!("{:12.4} ", his.block),
            "a/c" => match a_on_c {
                Some(a_on_c) => print!("{:12.4} ", a_on_c),
                None => print!("{:>12} ", "N/A"),
            },
            "a/d" => match his.geometry.a.as_ref() {
                Some(a) => print!("{:12.4} ", a.length / component.forward),
                None => print!("{:>12} ", "N/A"),
            },
            "c/b" => match his.geometry.c.as_ref() {
                Some(c) => print!("{:12.4} ", c.length / component.sideways),
                None => print!("{:>12} ", "N/A"),
            },
            "kmax_a" => match his.kmax.a {
                Some(kmax_a) => print!("{:12.4} ", kmax_a),
                None => print!("{:>12} ", "N/A"),
            },
            "kmax_c" => match his.kmax.c {
                Some(kmax_c) => print!("{:12.4} ", kmax_c),
                None => print!("{:>12} ", "N/A"),
            },
            "kmin_a" => match his.kmin.a {
                Some(kmin_a) => print!("{:12.4} ", kmin_a),
                None => print!("{:>12} ", "N/A"),
            },
            "kmin_c" => match his.kmin.c {
                Some(kmin_c) => print!("{:12.4} ", kmin_c),
                None => print!("{:>12} ", "N/A"),
            },
            "dk_a" => match his.dk.a {
                Some(dk_a) => print!("{:12.4} ", dk_a),
                None => print!("{:>12} ", "N/A"),
            },
            "dk_c" => match his.dk.c {
                Some(dk_c) => print!("{:12.4} ", dk_c),
                None => print!("{:>12} ", "N/A"),
            },
            "R" => if his.cycle.max.index == usize::MAX {
                print!("{:>12} ", "N/A")
            } else {
                print!("{:12.4} ", his.cycle.min.value / his.cycle.max.value)
            },
            "beta_a" => match his.beta.a {
                Some(beta_a) => print!("{:12.4e} ", beta_a),
                None => print!("{:>12} ", "N/A"),
            },
            "beta_c" => match his.beta.c {
                Some(beta_c) => print!("{:12.4e} ", beta_c),
                None => print!("{:>12} ", "N/A"),
            },
            "a" => match his.geometry.a.as_ref() {
                // Some(a) => print!("{} ", a.length),
                Some(a) => print!("{:12.6e} ", a.length),
                None => print!("{:>12} ", "N/A"),
            },
            "c" => match his.geometry.c.as_ref() {
                // Some(c) => print!("{} ", c.length),
                Some(c) => print!("{:12.6e} ", c.length),
                None => print!("{:>12} ", "N/A"),
            },
            "a_prior" => match his.geometry.a.as_ref() {
                Some(a) => print!("{:12.6e} ", a.length - his.length_deltas.a.unwrap_or(0.0)),
                None => print!("{:>12} ", "N/A"),
            },
            "c_prior" => match his.geometry.c.as_ref() {
                Some(c) => print!("{:12.6e} ", c.length - his.length_deltas.c.unwrap_or(0.0)),
                None => print!("{:>12} ", "N/A"),
            },
            "da" => match his.length_deltas.a {
                Some(da) => print!("{:12.4e} ", da),
                None => print!("{:>12} ", "N/A"),
            },
            "dc" => match his.length_deltas.c {
                Some(dc) => print!("{:12.4e} ", dc),
                None => print!("{:>12} ", "N/A"),
            },
            "peak" => if his.cycle.max.index == usize::MAX {
                print!("{:>12} ", "N/A")
            } else {
                print!("{:12.6} ", his.stress * his.cycle.max.value)
            },
            "valley" =>if his.cycle.max.index == usize::MAX {
                print!("{:>12} ", "N/A")
            } else {
                 print!("{:12.6} ", his.stress * his.cycle.min.value)
            },
            ref opt => {
                println!(
                    "Error: Unknown output option (use the --list option for a complete list): {}",
                    opt
                );
                process::exit(1);
            }
        };
    }
    println!();
}

pub fn format_parameters(
    method: Option<BTreeMap<OptimisableParameter, f64>>,
    dadn: BTreeMap<ParameterLabel, f64>
) -> String {
    let method = if method.is_some() {
        format!("{:#?}", method.unwrap())
    } else {
        "N/A".to_string()
    };

    format!("Method: {}\nda/dN: {:#?}", method, dadn)
}

/// This structure provides a way of remembering the failure messages
/// so that they can be written out at the end.
pub struct FailureResult {
    /// Type of failure.
    pub failure: bool,
    /// Failure Message.
    pub messages: String,
}

/// check if we have reached a pre-defined crack growth limit.
pub fn reached_limit(
    part_block: f64,
    block_limit: f64,
    geometry: &CrackGeometry,
    limit: &CrackLength,
) -> FailureResult {
    let mut message = String::new();
    let mut terminate = false;

    if part_block >= block_limit {
        let _ = writeln!(
            message,
            "{}Run stopped because hit block limit {}",
            COMMENT, block_limit
        );
        terminate = true;
    }

    if geometry.a.is_some() && limit.a.is_some() {
        let a = geometry.a.as_ref().unwrap().length;
        let a_limit = limit.a.unwrap();
        if a >= a_limit {
            let _ = writeln!(
                message,
                "{}Failure Event: a {:?} >= a_limit {:?}",
                COMMENT, a, a_limit
            );
            terminate = true;
        }
    }

    if geometry.c.is_some() && limit.c.is_some() {
        let c = geometry.c.as_ref().unwrap().length;
        let c_limit = limit.c.unwrap();
        if c >= c_limit {
            let _ = writeln!(
                message,
                "{}Failure Event: c {:?} >= c_limit {:?}",
                COMMENT, c, c_limit
            );
            terminate = true;
        }
    }

    FailureResult {
        failure: terminate,
        messages: message,
    }
}

/// Check whether the crack has exceeded any failure criteria for the component.
pub fn component_failed(
    geometry: &CrackGeometry,
    smax: f64,
    kmax: f64,
    component: &Component,
    optional_tests: &TerminationTests,
) -> FailureResult {
    let mut terminate = false;
    let mut message = "".to_string();

    // clippy says to do it this way, but its a bit inconsistent with the multiple failure checks
    // Check whether we have satisfied any termination criteria.
    // let mut terminate = if component.forward > 0.0 && a[0] > component.forward {
    //     message += &format!(
    //         "{}Failure Event: a[{}] > depth[{}]\n",
    //         COMMENT, a[0], component.forward
    //     );
    //     true } else { false };

    if geometry.a.is_some() {
        let a = geometry.a.as_ref().unwrap().length;
        if component.forward > 0.0 && a > component.forward {
            let _ = writeln!(
                message,
                "{}Failure Event: a[{}] > depth[{}]",
                COMMENT, a, component.forward
            );
            terminate = true;
        }
    }

    if geometry.c.is_some() {
        let c = geometry.c.as_ref().unwrap().length;
        if component.sideways > 0.0 && c > component.sideways {
            let _ = writeln!(
                message,
                "{}Failure Event: c[{}] > width[{}]",
                COMMENT, c, component.sideways
            );
            terminate = true;
        }
    }

    if optional_tests.limit_k1c && component.material.k1c > 0.0 && kmax > component.material.k1c {
        let _ = writeln!(
            message,
            "{}Failure Event: k[{}] > k1c[{}]",
            COMMENT, kmax, component.material.k1c
        );
        terminate = true;
    }

    if optional_tests.limit_yield {
        // The net stress for the component will depend on the shape of a crack.
        // We assume it is a corner crack (but the worst case will be an internal crack).
        let a = match &geometry.a {
            Some(a) => a.length,
            None => geometry.c.as_ref().unwrap().length,
        };

        let c = match &geometry.c {
            Some(c) => c.length,
            None => a,
        };

        let approx_crack_area = (PI / 4.0) * a * c;
        let component_area = component.sideways * component.forward;
        let applied_stress = smax * component_area / (component_area - approx_crack_area);

        if component.material.yield_stress > 0.0 && applied_stress > component.material.yield_stress
        {
            let _ = writeln!(
                message,
                "{}Note: Assuming a corner crack to check the net-section yield stress criterion.",
                COMMENT
            );
            let _ = writeln!(
                message,
                "{}approx crack area {}, component area {}",
                COMMENT, approx_crack_area, component_area
            );
            let _ = writeln!(
                message,
                "{}Failure Event: stress[{}] > yield[{}]",
                COMMENT, applied_stress, component.material.yield_stress
            );
            terminate = true;
        }
    }

    FailureResult {
        failure: terminate,
        messages: message,
    }
}
