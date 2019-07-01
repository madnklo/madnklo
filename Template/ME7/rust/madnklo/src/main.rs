extern crate libc;
extern crate madnklo;
extern crate rand;

use cuba::{CubaIntegrator, CubaVerbosity};

use madnklo::matrix_element_evaluator::{MatrixElementEvaluator, MatrixElement};
use madnklo::phase_space_generator::FlatPhaseSpaceGenerator;
use rand::Rng;
use libc::{c_char, c_double, c_int, c_longlong, c_void};

use madnklo::matrix_elements::LO_epem_ddxzz_1__1_epem_ddxzz_no_ememhzwpwpMatrixElement;

use std::f64::consts::PI;
use vector::LorentzVector;
struct Process {
    matrix_element: MatrixElementEvaluator<LO_epem_ddxzz_1__1_epem_ddxzz_no_ememhzwpwpMatrixElement>,
    phase_space_generator: FlatPhaseSpaceGenerator,
    external_momenta: Vec<LorentzVector<f64>>,
    e_cm: f64,
}

#[inline(always)]
fn integrand(
    x: &[f64],
    f: &mut [f64],
    process: &mut Process,
    nvec: usize,
    _core: i32,
) -> Result<(), &'static str> {
    let weight =
        process
            .phase_space_generator
            .generate(process.e_cm, x, &mut process.external_momenta[2..]);

    let res = process
        .matrix_element
        .evaluate(&process.external_momenta, -1, 0.118);

    if !res.is_finite() {
        println!("nan at x={:?}", x);
        for p in &process.external_momenta {
            println!("p={}, mass={:e}", p, p.square().abs().sqrt());
        }
    }

    f[0] = res * weight;

    // TODO: do after integration
    // multiply conversion factor to picobarns
    f[0] *= 0.389379304e9;

    // multiply the flux factor (with hardcoded 0 masses)
    f[0] /= (2. * process.e_cm * process.e_cm) * (2. * PI).powi(x.len() as i32);

    Ok(())
}

fn main() {
    let mut ci = CubaIntegrator::new(integrand);
    ci.set_nstart(100000)
        .set_nincrease(100000)
        .set_maxeval(10000000)
        .set_epsrel(0.0001)
        .set_seed(1)
        .set_cores(0, 1000);

    // hardcode e+ e- -> d d~ z z
    let e_cm = 1000.;
    let matrix_element =
        MatrixElementEvaluator::new("../Cards/param_card.dat", LO_epem_ddxzz_1__1_epem_ddxzz_no_ememhzwpwpMatrixElement {} );
    let fpsg = FlatPhaseSpaceGenerator::new(vec![0., 0., 91.188, 91.188]);
    let mut external_momenta = vec![LorentzVector::default(); 6]; // 2 incoming + 4 external
    external_momenta[0] = LorentzVector::from_args(e_cm / 2., 0., 0., e_cm / 2.);
    external_momenta[1] = LorentzVector::from_args(e_cm / 2., 0., 0., -e_cm / 2.);

    let data = Process {
        matrix_element,
        phase_space_generator: fpsg,
        external_momenta,
        e_cm,
    };
    let r = ci.vegas(8, 1, 1, CubaVerbosity::Progress, 0, data);

    println!("{:#?}", r);
}
