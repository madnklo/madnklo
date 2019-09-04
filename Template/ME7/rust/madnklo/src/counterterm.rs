use epsilon_expansion::EpsilonExpansion;
use std::collections::HashMap;
use std::f64::consts::PI;
use std::fmt;
use vector::LorentzVector;

pub trait CountertermEvaluator: Send + Sync {
    fn evaluate(
        &mut self,
        p: &HashMap<usize, LorentzVector<f64>>,
        alpha_s: f64,
        mu_r: f64,
        i_process: usize,
        i_CT: usize,
    ) -> vector<Event>;
}


pub struct Counterterm {

    // Ids
    process_id : i64,
    CT_id : i64,

    // Cards storing run, model and settings information
    param_card_path : str.
    run_card: RunCard,
    param_card: ParamCard,
    settings_card: SettingsCard,

    // Integrand evaluator containing all the hardcoded information
    counterterm_evaluator: Box<CountertermEvaluator>,

    // level of debug info
    verbosity: usize,
}

impl Counterterm {

    pub fn new(
        integrand_id : i64,
        process_id : i64,
        CT_id : i64,
        CT_mapping_id: i64,
        param_card_path: str,
        run_card: RunCard,
        param_card: ParamCard,
        settings_card: SettingsCard,
        counterterm_evaluator: Box<CountertermEvaluator>,
        subtraction_currents: Vec<SubtractionCurrent>,
    ) -> Integrand {

        Counterterm {
            integrand_id,
            process_id,
            CT_id,
            CT_mapping_id,
            param_card_path,
            run_card,
            param_card,
            settings_card,
            counterterm_evaluator,
            subtraction_currents: Vec<SubtractionCurrent>,
            verbosity: 0,
        }
    }

    pub fn set_verbosity(&mut self, verbosity: usize) {
        self.verbosity = verbosity;
    }

    pub fn evaluate(
        &mut self,
        process_id: usize,        
        PS_point: hashmap::<integer, LorenzVector>,
        base_weight: f64,
        mu_r: f64,
        mu_f1: f64,
        mu_f2: f64,
        xb_1: f64,
        xb_2: f64,
        x1: Option<f64>,
        x2: Option<f64>,
        all_resolved_flavors,
        hel_config : Option<i64>,
        apply_flavour_blind_cuts : bool,
        boost_back_to_com : bool,
        always_generate_event : bool,
        sector: Pair::<Option<Sector>,i64>,
    ) -> Vec<Event> {

        // STATIC NON-META-GENERATED CODE HERE
        let ME_call_ID = 1
        self.counterterm_evaluator.call_ME_for_key(
            &PS_point,
            self.param_card.aS,
            mu_r,
            process_id,
            CT_id,
            CT_mapping_id,
            ME_call_ID,
        )
        self.counterterm_evaluator.evaluate(
        )

        let SubtractionCurrentEvaluation = subtraction_currents[subtraction_current_call_ID].evaluate( ... )

        ...
        
    }

}

