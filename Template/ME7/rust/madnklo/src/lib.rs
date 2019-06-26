#[macro_use]
extern crate cpython;
#[macro_use]
extern crate lazy_static;
pub extern crate vector;

use cpython::PyResult;
use vector::LorentzVector;

pub mod phase_space_generator;

py_module_initializer!(madnklo, initmadnklo, PyInit_madnklo, |py, m| {
    m.add(py, "__doc__", "MadNkLO")?;
    m.add_class::<FlatPhaseSpaceGenerator>(py)?;
    Ok(())
});

py_class!(class FlatPhaseSpaceGenerator |py| {
    data gen: phase_space_generator::FlatPhaseSpaceGenerator;

    def __new__(_cls, masses: Vec<f64>) -> PyResult<FlatPhaseSpaceGenerator> {
        let gen = phase_space_generator::FlatPhaseSpaceGenerator::new(masses);
        FlatPhaseSpaceGenerator::create_instance(py, gen)
    }

    def generate(&self, e_cm: f64, r: Vec<f64>) -> PyResult<(Vec<Vec<f64>>, f64)> {
        // TODO: remove allocations
        let g = self.gen(py);
        let mut ps = vec![LorentzVector::default(); (r.len() + 4) / 3];
        let weight = g.generate(e_cm, &r, &mut ps);

        let psc = ps.iter().map(|p| vec![p.t, p.x, p.y, p.z]).collect();
        Ok((psc, weight))
    }
});