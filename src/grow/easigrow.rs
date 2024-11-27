//! The default easigrow growth method implementation

use std::{process, collections::BTreeMap};

use crate::{
    beta,
    cycle::Cycle,
    dadn::{self, ParameterLabel, relabel_parameters},
    tag::Tag,
};

use super::{Component, FatigueTest, History, BasicHistory, CrackLengthDelta, kmin, kmax, length_of_interest, Grow, CrackLength, StressIntensity, StressIntensityDelta};

#[derive(Clone)]
pub struct Easigrow {
    initial_history: History,
    history: Vec<History>,
    basic_history: Vec<BasicHistory>,
    grow: Grow,
    dadn: Box<dyn dadn::DaDn + Send + Sync>,
    beta: Box<dyn beta::Beta + Send + Sync>,
}

impl Easigrow {
    pub fn new(
        grow: Grow,
        initial: CrackLength,
        dadn: Box<dyn dadn::DaDn + Send + Sync>,
        beta: Box<dyn beta::Beta + Send + Sync>,
    ) -> Easigrow {
        let geometry = beta.geometry(initial.a.unwrap(), initial.c.unwrap());

        let beta_value = match beta.beta(&geometry) {
            Ok(beta_value) => beta_value,
            Err((_, message)) => {
                println!("Error: beta function failed with message '{}'. Check 
                starting crack lengths and component dimensions.", message);
                process::exit(1);
            }
        };

        let init_history = History::new(beta_value, geometry);

        Easigrow {
            initial_history: init_history.clone(),
            history: vec![init_history],
            basic_history: vec![],
            grow,
            dadn,
            beta,
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
            
            let kmax = match self.beta.direction_of_interest() {
                beta::DirectionOfInterest::A => moment.kmax.a.unwrap(),
                beta::DirectionOfInterest::C => moment.kmax.c.unwrap(),
            };

            let component_failed = super::component_failed(
                &moment.geometry,
                moment.stress * moment.cycle.max.value,
                kmax,
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

    fn grow_crack(&self, history: &History, cycle: &Cycle<Tag>, block: f64) -> History {
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

        let beta_result = match self.beta.beta(&history.geometry) {
            Ok(beta_result) => beta_result,
            Err((beta_result, _message)) => {
                beta_result
            }
        };

        // values around the crack front
        let mut length_deltas = CrackLengthDelta::default();
        let mut kmax_all = StressIntensity::default();
        let mut kmin_all = StressIntensity::default();
        let mut dk_all = StressIntensityDelta::default();
        let mut a_all = CrackLength::default();

        // Grow the crack in the 'a' direction, if necessary
        if let Some(beta_a) = beta_result.a {
            let a = history.geometry.a.as_ref().unwrap().length;
            let kmin = kmin(cycle, beta_a, a) * self.grow.scale;
            let kmax = kmax(cycle, beta_a, a) * self.grow.scale;
            let da = self.dadn.dadn(
                kmin,
                kmax,
                dadn::CrackState {
                    a,
                }
            );

            kmax_all.a = Some(kmax);
            kmin_all.a = Some(kmin);
            dk_all.a = Some(kmax - kmin);
            length_deltas.a = Some(da);
            a_all.a = Some(a + da);
            new_geometry.a.as_mut().unwrap().length = a_all.a.unwrap();
        }

        // Grow the crack in the 'c' direction, if necessary
        if let Some(beta_c) = beta_result.c {
            let c = history.geometry.c.as_ref().unwrap().length;
            let kmin = kmin(cycle, beta_c, c) * self.grow.scale;
            let kmax = kmax(cycle, beta_c, c) * self.grow.scale;
            let dc = self.dadn.dadn(
                kmin,
                kmax,
                dadn::CrackState {
                    a: c,
                }
            );

            kmax_all.c = Some(kmax);
            kmin_all.c = Some(kmin);
            dk_all.c = Some(kmax - kmin);
            length_deltas.c = Some(dc);
            a_all.c = Some(c + dc);
            new_geometry.c.as_mut().unwrap().length = a_all.c.unwrap();
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
        }
    }

    fn decode_parameter_slice(&self, parameters: &[f64]) -> BTreeMap<ParameterLabel, f64> {
        relabel_parameters(parameters, self.dadn.get_name()).unwrap()
    }

    fn reset(&mut self) {
        self.history = vec![self.initial_history.clone()];

        let crack_length = length_of_interest(self.beta.direction_of_interest(), &self.initial_history.geometry);

        let basic_history = BasicHistory {
            block: 0.0,
            crack_length,
        };
        self.basic_history = vec![basic_history];
    }

    fn printable_parameters(&self, parameters: &[f64]) -> String {
        let dadn_params = relabel_parameters(parameters, self.dadn.get_name()).unwrap();

        super::format_parameters(None, dadn_params)
    }
}

impl FatigueTest for Easigrow {
    fn run(&mut self) -> (&[History], String) {
        let message = self.run(true);

        (&self.history, message)
    }

    fn run_for_optimisation(&mut self) -> &[BasicHistory] {
        let _ = self.run(false);
        
        &self.basic_history
    }

    fn update_parameters(&mut self, parameters: &[f64]) {
        let dadn_params = self.decode_parameter_slice(parameters);
        self.dadn.update_parameters(&dadn_params, None);
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
