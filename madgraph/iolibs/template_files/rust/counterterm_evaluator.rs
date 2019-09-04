
pub struct CountertermEvaluator_%(counterterm_short_name)s {
    integrand_id; usize,
    process_id : usize,
    counterterm_id: usize,
    mapping_id: usize,
    %(counterterm_evaluator_definition)s
}

impl CountertermEvaluator_%(counterterm_short_name)s {

    pub fn new(card_filename: &str) -> CountertermEvaluator_%(counterterm_short_name)s {
        CountertermEvaluator_%(counterterm_short_name)s {
            integrand_id : %(integrand_id)d,
            process_id : %(process_id)d,
            counterterm_id: %(counterterm_id)d,
            mapping_id: %(mapping_id)d,
            %(counterterm_evaluator_construction)s,
        }
    }

}

impl CountertermEvaluator for CountertermEvaluator_%(counterterm_short_name)s {

    pub fn evaluate(
        &mut self,
        PS_point: hashmap::<integer, LorenzVector>,
        base_weight: f64,
        mu_r: f64,
        mu_f1: f64,
        mu_f2: f64,
        xb_1: f64,
        xb_2: f64,
        x1: Option<f64>,
        x2: Option<f64>,
        hel_config : Option<i64>,
        apply_flavour_blind_cuts : bool,
        boost_back_to_com : bool,
        always_generate_event : bool,
        sector: Pair::<Option<Sector>,i64>,
    ) -> Vec<Event> {

        // DYNAMICALLY META-GENERATED CODE HERE        
        %(counterterm_evaluation_code)s

    }

    fn call_ME_for_key(&mut self, p : &HashMap<usize , LorentzVector<f64>>, spin_correlation_vectors : &Vec<LorentzVector<f64>>, alpha_s : f64, mu_r : f64, i_process: usize, i_CT: usize, i_mapping: usize, i_call: usize) -> EpsilonExpansion<U13, N6> {
        match (i_process, i_CT, i_mapping, i_call) {
            %(ME_calls)s
            _ => panic!("The matrix element call signature ({},{},{},{}) for integrand id {} is not implemented.", i_process, i_CT, i_mapping, i_call, self.integrand_id)
        }
    }

}
