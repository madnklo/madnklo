%(header)s

struct IntegrandEvaluator_%(integrand_ID)d {
    ID : usize,
}

impl IntegrandEvaluator_%(integrand_ID)d {

    fn new() -> IntegrandEvaluator_%(integrand_ID)d {
        IntegrandEvaluator_%(integrand_ID)d {
            %(integrand_ID)d
        }
    }

    fn call_ME_for_key(&mut self, p : HashMap<usize : LorentzVector>, alpha_s : f64, mu_r : f64, i_process: usize, i_CT: usize, i_call: usize) {

        return match (i_process, i_CT, i_call) {
            %(ME_calls)s
            _ => panic!("The matrix element call signature ({},{},{}) for integrand ID {} is not implemented.".format!(i_process, i_CT, i_call, self.ID))
        }
    }

}
