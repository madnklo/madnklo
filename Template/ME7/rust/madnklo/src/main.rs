extern crate libc;
extern crate madnklo;
extern crate rand;

use cuba::{CubaIntegrator, CubaVerbosity};

use libc::{c_char, c_double, c_int, c_longlong, c_void};
use madnklo::integrand::Integrand;
use madnklo::matrix_element_evaluator::{MatrixElement, MatrixElementEvaluator};
use madnklo::phase_space_generator::{FlatPhaseSpaceGenerator, PhaseSpaceGenerator};
use rand::Rng;

use madnklo::all_integrands;

use std::f64::consts::PI;
use vector::LorentzVector;
struct Process<'a> {
    integrand: &'a mut Integrand,
}

#[inline(always)]
fn integrand(
    x: &[f64],
    f: &mut [f64],
    process: &mut Process,
    nvec: usize,
    _core: i32,
) -> Result<(), &'static str> {
    let res = process.integrand.evaluate(x, 1.0, None);

    if !res.is_finite() {
        println!("nan at x={:?}", x);
    }

    f[0] = res;

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

    let mut integrands = all_integrands::AllIntegrands_RUSTTEST::new();
    let integrand = integrands.all_integrands.get_mut(&1).unwrap();
    integrand.set_verbosity(1);

    // for testing
    println!(
        "Sample point #1: {}",
        integrand.evaluate(
            &[
                0.8072168371449964,
                0.09703446753764455,
                0.3807885551627419,
                0.15528816294348846,
                0.43059185758737006,
                0.09805328339852282,
                0.8104023752957474,
                0.7673972319392628,
            ],
            1.0,
            None,
        )
    );

    println!(
        "Sample point #2: {}",
        integrand.evaluate(
            &[
                0.5078178077029998,
                0.7861861126143361,
                0.1830110046211958,
                0.9629787226544925,
                0.003683355541845512,
                0.24911371478768107,
                0.6216002567286699,
                0.2602448561130817,
            ],
            1.0,
            None,
        )
    );

    integrand.set_verbosity(0);
    let dims = integrand.get_dimensions();

    let data = Process {
        integrand: integrand,
    };
    let r = ci.vegas(dims, 1, 1, CubaVerbosity::Progress, 0, data);

    println!("{:#?}", r);
}
