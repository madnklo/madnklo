#[macro_use]
extern crate cpython;
#[macro_use]
extern crate lazy_static;
extern crate cuba;
pub extern crate vector;
pub extern crate epsilon_expansion;

use cpython::PyResult;
use vector::LorentzVector;
use std::cell::RefCell;


macro_rules! hashmap {
    ($( $key: expr => $val: expr ),*) => {{
         let mut map = ::std::collections::HashMap::new();
         $( map.insert($key, $val); )*
         map
    }}
}

pub mod matrix_element_evaluator;
pub mod phase_space_generator;
pub mod matrix_elements; // generated matrix elements
pub mod integrand;
pub mod run_card;
pub mod param_card;
pub mod settings_card;
pub mod all_integrands;
pub mod integrands;

use crate::phase_space_generator::PhaseSpaceGenerator;

py_module_initializer!(madnklo, initmadnklo, PyInit_madnklo, |py, m| {
    m.add(py, "__doc__", "MadNkLO")?;
    m.add_class::<FlatPhaseSpaceGenerator>(py)?;
    Ok(())
});

py_class!(class FlatPhaseSpaceGenerator |py| {
    data gen: RefCell<phase_space_generator::FlatPhaseSpaceGenerator>;

    def __new__(_cls, masses: Vec<f64>) -> PyResult<FlatPhaseSpaceGenerator> {
        let gen = phase_space_generator::FlatPhaseSpaceGenerator::new(masses);
        FlatPhaseSpaceGenerator::create_instance(py, RefCell::new(gen))
    }

    def generate(&self, e_cm: f64, r: Vec<f64>) -> PyResult<(Vec<Vec<f64>>, f64)> {
        // TODO: remove allocations
        let mut g = self.gen(py).borrow_mut();
        let mut ps = vec![LorentzVector::default(); (r.len() + 4) / 3];
        let weight = g.generate(e_cm, &r, &mut ps);

        let psc = ps.iter().map(|p| vec![p.t, p.x, p.y, p.z]).collect();
        Ok((psc, weight))
    }
});