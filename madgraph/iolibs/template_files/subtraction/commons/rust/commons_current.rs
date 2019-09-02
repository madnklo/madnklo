use epsilon_expansion::EpsilonExpansion;
use std::collections::HashMap;
use std::f64::consts::PI;
use std::fmt;
use vector::LorentzVector;

pub trait SubtractionCurrentEvaluator: Send + Sync {
    fn evaluate(
        &mut self,
        p: &HashMap<usize, LorentzVector<f64>>,
        alpha_s: f64,
        mu_r: f64,
        i_process: usize,
        i_CT: usize,
        i_current: usize,
    ) -> vector<Event>;
}


pub struct SubtractionCurrent {

    // Cards storing run, model and settings information
    run_card: RunCard,
    param_card: ParamCard,
    settings_card: SettingsCard,

    // Integrand evaluator containing all the hardcoded information
    subtraction_current_evaluator: Box<SubtractionCurrentEvaluator>,

    // level of debug info
    verbosity: usize,
}

impl SubtractionCurrent {

    pub fn new(
        param_card: ParamCard,
        subtraction_settings_card: SettingsCard,
        subtraction_current_evaluator: Box<SubtractionCurrentEvaluator>,
    ) -> SubtractionCurrent {

        SubtractionCurrent {
            run_card,
            param_card,
            settings_card,
            subtraction_current_evaluator,
            verbosity: 0,
        }
    }

    pub fn set_verbosity(&mut self, verbosity: usize) {
        self.verbosity = verbosity;
    }

}

