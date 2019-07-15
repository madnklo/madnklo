use std::collections::HashMap;
use vector::LorentzVector;
use epsilon_expansion::EpsilonExpansion;
use typenum::{U13,N6};
use crate::integrand::IntegrandEvaluator;

use crate::matrix_element_evaluator::MatrixElementEvaluator;
%(header)s

pub struct IntegrandEvaluator_%(integrand_short_name)s {
    ID : usize,
    %(integrand_evaluator_definition)s
}

impl IntegrandEvaluator_%(integrand_short_name)s {

    pub fn new(card_filename: &str) -> IntegrandEvaluator_%(integrand_short_name)s {
        IntegrandEvaluator_%(integrand_short_name)s {
            ID: %(integrand_ID)d,
            %(integrand_evaluator_construction)s
        }
    }

}

impl IntegrandEvaluator for IntegrandEvaluator_%(integrand_short_name)s {

    fn call_ME_for_key(&mut self, p : &HashMap<usize , LorentzVector<f64>>, alpha_s : f64, mu_r : f64, i_process: usize, i_CT: usize, i_call: usize) -> EpsilonExpansion<U13, N6> {
        match (i_process, i_CT, i_call) {
            %(ME_calls)s
            _ => panic!("The matrix element call signature ({},{},{}) for integrand ID {} is not implemented.", i_process, i_CT, i_call, self.ID)
        }
    }

}
