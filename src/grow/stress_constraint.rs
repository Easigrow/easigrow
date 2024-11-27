use std::collections::HashMap;

#[derive(Debug, Clone, Copy)]
pub struct StressConstraints {
    pub variable_alpha: Option<VariableAlpha>,
    pub alpha: Option<f64>,
}

pub enum StressConstraint {
    PlasticTension,
    PlasticCompression,
    WakeCompression,
}

#[derive(Debug, Clone, Copy)]
pub enum CalculationMode {
    GrowthRate,
    LogGrowthRate,
    CrackLength,
    LogCrackLength,
}

pub struct VariableAlphaInputs {
    pub growth_rate: f64,
    pub crack_length: f64,
}

#[derive(Debug, Clone, Copy)]
pub struct VariableAlpha {
    pub alpha_max: f64,
    pub alpha_min: f64,
    // Transition rate
    pub k: f64,
    // Centre point
    pub x0: f64,
    pub calculation_mode: CalculationMode,
}

impl VariableAlpha {
    pub fn from_map(
        data: &HashMap<String, f64>,
        calculation_mode: CalculationMode,
    ) -> Result<Self, String> {
        if data.len() != 4 {
            return Err("Incorrect number of parameters given; expecting 4.".to_string());
        }

        let mut alpha_max = 0.0;
        let mut alpha_min = 0.0;
        let mut k = 0.0;
        let mut x0 = 0.0;

        for (key, value) in data {
            match key.as_str() {
                "alpha_max" => alpha_max = *value,
                "alpha_min" => alpha_min = *value,
                "k" => k = *value,
                "x0" => x0 = *value,
                _ => return Err("Invalid variable constraint parameter: ".to_string() + key),
            }
        }

        Ok(Self {
            alpha_max,
            alpha_min,
            k,
            x0,
            calculation_mode,
        })
    }

    pub fn calculate(&self, inputs: &VariableAlphaInputs) -> f64 {
        let input_value = match self.calculation_mode {
            CalculationMode::GrowthRate => inputs.growth_rate,
            CalculationMode::LogGrowthRate => inputs.growth_rate.ln(),
            CalculationMode::CrackLength => inputs.crack_length,
            CalculationMode::LogCrackLength => inputs.crack_length.ln(),
        };

        self.alpha_max
            - (self.alpha_max - self.alpha_min) / (1.0 + (-self.k * (input_value - self.x0)).exp())
    }
}
