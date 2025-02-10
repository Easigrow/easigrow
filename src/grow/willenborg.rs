//! Implementation of the Willenborg crack growth method

use std::{collections::BTreeMap, process};

use crate::{
    beta,
    cycle::Cycle,
    dadn::{self, ParameterLabel, relabel_parameters},
    plastic::zone_size_new,
    tag::Tag,
};

use super::{kmax, kmin, length_of_interest, stress_constraint::{VariableAlpha, VariableAlphaInputs}, BasicHistory, Component, CrackLength, CrackLengthDelta, FatigueTest, Grow, History, OptimisableParameter, StressIntensity, StressIntensityDelta};

#[derive(Clone)]
struct Overload {
    dmax: f64,
    dist_plz_ol: f64,
    kmax_ol: f64,
}

#[derive(Clone, Default)]
struct Overloads {
    a: Overload,
    c: Overload,
}

#[derive(Clone)]
pub struct Willenborg {
    initial_history: History,
    history: Vec<History>,
    basic_history: Vec<BasicHistory>,
    grow: Grow,
    dadn: Box<dyn dadn::DaDn + Send + Sync>,
    beta: Box<dyn beta::Beta + Send + Sync>,
    alpha: Option<f64>,
    variable_alpha: Option<VariableAlpha>,
    overloads: Overloads,
    deltak_th: f64,
    r_so: f64,
}

impl Willenborg {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        grow: Grow,
        initial: CrackLength,
        dadn: Box<dyn dadn::DaDn + Send + Sync>,
        beta: Box<dyn beta::Beta + Send + Sync>,
        alpha: Option<f64>,
        variable_alpha: Option<VariableAlpha>,
        deltak_th: f64,
        r_so: f64,
    ) -> Willenborg {
        let geometry = beta.geometry(initial.a.unwrap(), initial.c.unwrap());

        let beta_result = match beta.beta(&geometry) {
            Ok(beta_result) => beta_result,
            Err((_, message)) => {
                println!("Error: beta function failed with message '{}'. Check 
                starting crack lengths and component dimensions.", message);
                process::exit(1);
            }
        };

        let init_history = History::new(beta_result, geometry);

        Willenborg {
            initial_history: init_history.clone(),
            history: vec![init_history],
            basic_history: vec![],
            grow,
            dadn,
            beta,
            alpha,
            variable_alpha,
            overloads: Overloads::default(),
            deltak_th,
            r_so,
        }
    }

    fn run(&mut self, save_full_history: bool) -> String {
        self.reset();
        let message;
        let mut block = 0.0;
        let mut position = 0;
        let mut previous_history = self.initial_history.clone();
        let mut keep_every_counter = 0;
        
        loop {
            let cycle = self.grow.cycles[position];
            let moment = self.grow_crack(&previous_history, &cycle, block);

            // Check for terminating conditions.
            let reached_limit = super::reached_limit(
                block,
                self.grow.block_limit,
                &moment.geometry,
                &self.grow.limit,
            );
            
            let component_failed = super::component_failed(
                &moment.geometry,
                moment.stress * moment.cycle.max.value,
                moment.kmax.c.unwrap(),
                &self.grow.component,
                &self.grow.termination_tests,
            );

            keep_every_counter += 1;

            if save_full_history {
                self.history.push(moment.clone());
            } else if keep_every_counter == self.grow.keep_every {
                let crack_length = length_of_interest(self.beta.direction_of_interest(), &moment.geometry);

                let basic_moment = BasicHistory {
                    block,
                    crack_length,
                };

                self.basic_history.push(basic_moment);
                keep_every_counter = 0;
            }
            
            if reached_limit.failure || component_failed.failure {
                message = reached_limit.messages + &component_failed.messages;
                break;
            }

            block = block.floor() + (position as f64 + 1.0) / self.grow.cycles.len() as f64;
            previous_history = moment;
            position = (position + 1) % self.grow.cycles.len();
        }

        message
    }

    fn grow_crack(&mut self, history: &History, cycle: &Cycle<Tag>, block: f64) -> History {
        // let history = self.history.last().unwrap();
        let scale = self.grow.scale;
        let smax = cycle.max.value;
        let smin = cycle.min.value;
        let mut new_geometry = history.geometry.clone();

        if smax < smin {
            println!(
                "Program Error: smax {} is less than smin {}. This should never happen.",
                smax, smin
            );
            process::exit(1);
        }

        // TODO: Properly handle failing beta
        let beta_result = match self.beta.beta(&history.geometry) {
            Ok(beta_result) => beta_result,
            Err((beta_result, _message)) => {
                beta_result
            }
        };

        // Values around the crack front
        let mut length_deltas = CrackLengthDelta::default();
        let mut kmax_all = StressIntensity::default();
        let mut kmin_all = StressIntensity::default();
        let mut dk_all = StressIntensityDelta::default();
        let mut a_all = CrackLength::default();

        let mut kmax_a = 0.0;
        let mut kmin_a = 0.0;
        let mut kmax_c = 0.0;
        let mut kmin_c = 0.0;
        let mut a = 0.0;
        let mut c = 0.0;

        if let Some(beta_c) = beta_result.c {
            c = match history.geometry.get_length_c() {
                Some(length) => length,
                None => {
                    println!("Error: Expected length 'c' in geometry, none found.");
                    process::exit(1);
                }
            };
            kmax_c = kmax(cycle, beta_c, c) * self.grow.scale;
            kmin_c = kmin(cycle, beta_c, c) * self.grow.scale;
        }

        if let Some(beta_a) = beta_result.a {
            a = match history.geometry.get_length_a() {
                Some(length) => length,
                None => {
                    println!("Error: Expected length 'a' in geometry, none found.");
                    process::exit(1);
                }
            };
            kmax_a = kmax(cycle, beta_a, a) * self.grow.scale;
            kmin_a = kmin(cycle, beta_a, a) * self.grow.scale;
        }

        let alpha = match self.beta.direction_of_interest() {
            beta::DirectionOfInterest::A => {
                let growth_rate = self.dadn.dadn(kmin_a, kmax_a, dadn::CrackState{a});
                self.alpha(a, growth_rate)
            },
            beta::DirectionOfInterest::C => {
                let growth_rate = self.dadn.dadn(kmin_c, kmax_c, dadn::CrackState{a: c});
                self.alpha(c, growth_rate)
            },
        };
        
        if beta_result.c.is_some() {
            let (delta_c, new_overload) = self.delta(kmin_c, kmax_c, c, &self.overloads.c, alpha);
            let new_length = c + delta_c;
            length_deltas.c = Some(delta_c);
            a_all.c = Some(new_length);
            kmax_all.c = Some(kmax_c);
            kmin_all.c = Some(kmin_c);
            dk_all.c = Some(kmax_c - kmin_c);
            new_geometry.set_length_c(new_length);
    
            if let Some(overload) = new_overload {
                self.overloads.c = overload;
            }
        }

        if beta_result.a.is_some() {
            let (delta_a, new_overload) = self.delta(kmin_a, kmax_a, a, &self.overloads.a, alpha);
            let new_length = a + delta_a;
            length_deltas.a = Some(delta_a);
            a_all.a = Some(new_length);
            kmax_all.a = Some(kmax_a);
            kmin_all.a = Some(kmin_a);
            dk_all.a = Some(kmax_a - kmin_a);
            new_geometry.set_length_a(new_length);

            if let Some(overload) = new_overload {
                self.overloads.a = overload;
            }
        }

        History {
            block, // part_block,
            length_deltas,
            kmax: kmax_all,
            kmin: kmin_all,
            dk: dk_all,
            cycle: *cycle,
            stress: scale,
            beta: beta_result,
            geometry: new_geometry,
            r: smin / smax,
            peak: scale * cycle.max.value,
            valley: scale * cycle.min.value,
        }
    }

    /// Calculate the delta for a crack front.
    /// 
    /// Returns a tuple of `(f64, Option<Overload>)`, containing the result and a new
    /// overload if one is produced.
    fn delta(&self, mut kmin: f64, mut kmax: f64, crack_length: f64, overload: &Overload, alpha: f64) -> (f64, Option<Overload>) {
        if kmax < self.deltak_th {
            // When kmax < deltak_th, the result is erroneous phi values.
            // So don't calculate any plastic zone or growth in this case.
            (0.0, None)
        } else {
            // Calculate the size of the plastic zone
            let plastic_zone = zone_size_new(
                kmax,
                self.grow.component.material.yield_stress,
                alpha,
            );
            let d = crack_length + plastic_zone;
            let mut new_overload = None;

            if d >= overload.dmax {
                // New overload, calculate without any modifications
                new_overload = Some(Overload {
                    dmax: d,
                    dist_plz_ol: plastic_zone,
                    kmax_ol: kmax,
                });
            } else {
                // Calculate effective kmax and kmin
                let phi = (1.0 - (self.deltak_th / kmax)) / (self.r_so - 1.0);
                let k_ap = overload.kmax_ol * ((overload.dmax - crack_length) / overload.dist_plz_ol).sqrt();
                let k_red = phi * (kmax - k_ap);
                kmin += k_red;
                kmax += k_red;
            }

            let delta = if kmax >= 0.0 {
                self.dadn.dadn(
                    kmin,
                    kmax,
                    dadn::CrackState {
                        a: crack_length,
                    }
                )
            } else {
                0.0
            };

            (delta, new_overload)
        }
    }

    fn alpha(&self, crack_length: f64, growth_rate: f64) -> f64 {
        if let Some(constant_value) = self.alpha {
            constant_value
        } else {
            let variable_alpha_inputs = VariableAlphaInputs {
                crack_length,
                growth_rate,
            };
            self.variable_alpha.unwrap().calculate(&variable_alpha_inputs)
        } 
    }

    // Ensure that the order of the parameters defined here matches the order in 
    // decode_parameter_slice().
    pub fn optimisable_parameters(is_variable_alpha: bool) -> Vec<OptimisableParameter> {
        use OptimisableParameter::*;

        if is_variable_alpha{
            vec![rso, deltak_th, alpha_min, alpha_max, k, x0]
        } else {
            vec![rso, deltak_th, alpha]
        }
    }

    // Ensure that the order of the parameters defined here matches the order in 
    // optimisable_parameters().
    fn decode_parameter_slice(&mut self, parameters: &[f64]) {
        let dadn_start;
        self.r_so = parameters[0];
        self.deltak_th = parameters[1]; // This may be overwritten below

        if let Some(ref mut variable_alpha) = self.variable_alpha {
            variable_alpha.alpha_min = parameters[2];
            variable_alpha.alpha_max = parameters[3];
            variable_alpha.k = parameters[4];
            variable_alpha.x0 = parameters[5];
            dadn_start = 6;
        } else {
            self.alpha = Some(parameters[2]);
            dadn_start = 3;
        }
        
        let dadn_params = relabel_parameters(&parameters[dadn_start..], self.dadn.get_name()).unwrap();

        // If a deltak_th is used in the dadn, it should be the same as the deltak_th used in the method.
        // Use the dadn's deltak_th given by the optimiser and ignore the one allocated to the method.
        if let Some(deltak_th) = dadn_params.get(&ParameterLabel::deltak_th) {
            self.deltak_th = *deltak_th;
        }
        
        self.dadn.update_parameters(&dadn_params, Some(self.deltak_th));
    }

    fn printable_parameters(&self, parameters: &[f64]) -> String {
        let mut method_params = BTreeMap::new();

        for (i, param) in Willenborg::optimisable_parameters(self.variable_alpha.is_some()).iter().enumerate() {
            method_params.insert(*param, parameters[i]);
        }

        let dadn_start = method_params.len();
        let dadn_params = relabel_parameters(&parameters[dadn_start..], self.dadn.get_name()).unwrap();

        // The parameter list contains an unused deltak_th for the method in the instances where
        // a dadn includes a deltak_th. The data handling for deltak_th happens in decode_parameter_slice
        if let Some(deltak_th) = dadn_params.get(&ParameterLabel::deltak_th) {
            method_params.insert(OptimisableParameter::deltak_th, *deltak_th);
        }

        super::format_parameters(Some(method_params), dadn_params)
    }

    fn reset(&mut self) {
        self.history = vec![self.initial_history.clone()];

        let crack_length = length_of_interest(self.beta.direction_of_interest(), &self.initial_history.geometry);

        let basic_history = BasicHistory {
            block: 0.0,
            crack_length,
        };
        self.basic_history = vec![basic_history];

        self.overloads = Overloads::default();
    }
}

impl FatigueTest for Willenborg {
    fn run(&mut self) -> (&[History], String) {
        let message = self.run(true);

        (&self.history, message)
    }

    fn run_for_optimisation(&mut self) -> &[BasicHistory] {
        let _ = self.run(false);
        
        &self.basic_history
    }

    fn update_parameters(&mut self, parameters: &[f64]) {
        self.decode_parameter_slice(parameters);
    }

    fn reset(&mut self) {
        self.reset();
    }

    fn get_component(&self) -> &Component {
        &self.grow.component
    }

    fn get_cycles(&self) -> &[Cycle<Tag>] {
        &self.grow.cycles
    }

    fn get_history(&self) -> &[History] {
        &self.history
    }

    fn inner_clone(&self) -> Box<dyn FatigueTest + Send + Sync> {
        Box::new(self.clone())
    }

    fn printable_parameters(&self, parameters: &[f64]) -> String {
        self.printable_parameters(parameters)
    }
}

// pub fn initialise_cycles(sequence: &[Tag]) -> (Vec<Cycle<Tag>>, Vec<Tag>) {
//     let mut turning_points = cycle::turning_points(sequence);
//     cycle::balance_ends(&mut turning_points);
//     cycle::start_with_valley(&mut turning_points);

//     // Rainflow count
//     let (mut cycles, leftover) = cycle::unsorted_rainflow(&turning_points);

//     // Tension count the leftovers and add to the rest
//     let mut tension;
//     let mut unclosed = vec![];
//     if !leftover.is_empty() {
//         (tension, unclosed) = cycle::tension(&leftover);
//         cycles.append(&mut tension);
//     }

//     cycle::sort_peak_order(&mut cycles);

//     (cycles, unclosed)
// }

impl Default for Overload {
    fn default() -> Self {
        Self {
            dmax: 0.0,
            dist_plz_ol: 0.0,
            kmax_ol: 0.0
        }
    }
}

#[cfg(test)]
mod tests {
    use std::f64::consts::FRAC_PI_2;

    use crate::{beta, cycle::Cycle, dadn, grow::{willenborg::Overloads, Component, CrackFront, CrackGeometry, CrackLength, Grow, History, TerminationTests}, material, tag::Tag};

    use super::Willenborg;

    fn three_simple_cycles() -> Vec<Cycle<Tag>> {
        let mut cycles = Vec::new();
        let sequence = vec![0.0, 8.9, 0.0, 5.4];
        let sequence = Tag::from(&sequence);

        for i in (0..sequence.len()).step_by(2) {
            cycles.push(Cycle{min: sequence[i], max: sequence[i + 1]});
        }

        cycles
    }

    #[test]
    fn run_returns_correct_history_for_three_cycle_input() {
        let component = Component {
            sideways: 1.0,
            forward: f64::INFINITY,
            radius: f64::INFINITY,
            material: material::Properties::default(),
        };

        let termination_tests = TerminationTests {
            limit_k1c: true,
            limit_yield: true,
        };

        let cycles = three_simple_cycles();

        let a = 0.0;
        let c = 1e-3;
        let initial = CrackLength{ a: Some(a), c: Some(c) };
        let limit = CrackLength{ a: Some(1.0), c: Some(1.0) };
        let geometry = CrackGeometry {
            a: Some(CrackFront {
                length: a,
                angle: 0.0,
            }),
            c: Some(CrackFront {
                length: c,
                angle: FRAC_PI_2,
            }),
            ratio: Some(a / c),
        };
        let beta = beta::get_beta_fn("ct-fedderson66", &component, geometry, None);

        let dadn_params = &material::get_dadn("paris:default").unwrap().params;
        let dadn_options = dadn::Options{deltak_th: 0.0, rmax: 1.0, rmin: -1.0, kneg: false, kmax_th: 0.0, kmax_th_follows_deltak_th: true };
        let dadn = dadn::make_model("paris", dadn_params, dadn_options).unwrap();

        let geometry = beta.geometry(initial.a.unwrap(), initial.c.unwrap());
        let beta_result = beta.beta(&geometry).unwrap();

        let init_history = History::new(beta_result, geometry);

        let grow = Grow {
            cycles,
            component,
            scale: 10.0,
            limit,
            block_limit: 1.0,
            termination_tests,
            keep_every: 1,
        };

        let mut method = Willenborg {
            initial_history: init_history.clone(),
            history: vec![init_history],
            basic_history: vec![],
            grow,
            dadn,
            beta,
            alpha: Some(1.0),
            variable_alpha: None,
            overloads: Overloads::default(),
            deltak_th: 0.69,
            r_so: 3.0,
        };

        let message = method.run(true);

        println!("{:#?}", method.history);
        println!("{}", message);

        ()
    }
}
