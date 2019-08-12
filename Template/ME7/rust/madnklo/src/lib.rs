#![allow(non_snake_case)]
#![allow(non_camel_case_types)]

#[macro_use]
extern crate cpython;
extern crate cuba;
pub extern crate epsilon_expansion;
pub extern crate vector;

use cpython::PyResult;
use std::cell::RefCell;
use std::hash;
use std::mem;
use vector::LorentzVector;

macro_rules! hashmap {
    ($( $key: expr => $val: expr ),*) => {{
         let mut map = ::std::collections::HashMap::new();
         $( map.insert($key, $val); )*
         map
    }}
}

#[derive(Debug, Copy, Clone)]
struct HashableFloat(f64);

impl HashableFloat {
    fn key(&self) -> u64 {
        unsafe { mem::transmute(self.0) }
    }
}

impl hash::Hash for HashableFloat {
    fn hash<H>(&self, state: &mut H)
    where
        H: hash::Hasher,
    {
        self.key().hash(state)
    }
}

impl PartialEq for HashableFloat {
    fn eq(&self, other: &HashableFloat) -> bool {
        self.key() == other.key()
    }
}

impl Eq for HashableFloat {}

pub mod all_integrands;
pub mod integrand;
pub mod integrands;
pub mod matrix_element_evaluator;
pub mod matrix_elements; // generated matrix elements
pub mod param_card;
pub mod phase_space_generator;
pub mod run_card;
pub mod settings_card;

use crate::phase_space_generator::PhaseSpaceGenerator;

py_module_initializer!(madnklo, initmadnklo, PyInit_madnklo, |py, m| {
    m.add(py, "__doc__", "MadNkLO")?;
    m.add_class::<FlatPhaseSpaceGenerator>(py)?;
    Ok(())
});

py_class!(class FlatPhaseSpaceGenerator |py| {
    data gen: RefCell<phase_space_generator::FlatPhaseSpaceGenerator>;

    def __new__(_cls, n_initial: usize,
            masses: (Vec<f64>, Vec<f64>),
            collider_energy: f64,
            beam_type: (isize, isize),
            correlated_beam_convolution: bool,
            is_beam_factorization_active: (bool, bool)) -> PyResult<FlatPhaseSpaceGenerator> {
        let gen = phase_space_generator::FlatPhaseSpaceGenerator::new(n_initial, masses,
                collider_energy, beam_type, correlated_beam_convolution, is_beam_factorization_active);
        FlatPhaseSpaceGenerator::create_instance(py, RefCell::new(gen))
    }

    def generate(&self, e_cm: f64, r: Vec<f64>) -> PyResult<(Vec<Vec<f64>>, f64)> {
        // TODO: remove allocations
        let g = self.gen(py).borrow_mut();
        let mut ps = vec![LorentzVector::default(); (r.len() + 4) / 3];
        let weight = g.generate(e_cm, &r, &mut ps);

        let psc = ps.iter().map(|p| vec![p.t, p.x, p.y, p.z]).collect();
        Ok((psc, weight))
    }

    def get_PS_point(&self, x: Vec<f64>) -> PyResult<(Vec<Vec<f64>>, f64, (f64, Option<f64>), (f64, Option<f64>))> {
        // TODO: remove allocations
        let mut g = self.gen(py).borrow_mut();
        let mut ps = vec![LorentzVector::default(); g.masses.0.len() + g.masses.1.len()];

        let (wgt, (xb_1, xi1), (xb_2, xi2)) = g.get_PS_point(&x, &mut ps);

        let psc = ps.iter().map(|p| vec![p.t, p.x, p.y, p.z]).collect();
        Ok((psc, wgt, (xb_1, xi1), (xb_2, xi2)))
    }

});