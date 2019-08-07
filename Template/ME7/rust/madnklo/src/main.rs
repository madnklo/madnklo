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
    let dims = integrand.get_dimensions();

    integrand.set_verbosity(1);

    // for testing
    let point1 = [
        0.8072168371449964,
        0.09703446753764455,
        0.3807885551627419,
        0.15528816294348846,
        0.43059185758737006,
        0.09805328339852282,
        0.8104023752957474,
        0.7673972319392628,
        0.3196893458725276,
        0.7613783440709990,
        0.4689202598577768,
        0.1085951805298387,
        0.5348251805222444,
    ];

    let point2 = [
        0.8072168371449964,
        0.09703446753764455,
        0.3807885551627419,
        0.15528816294348846,
        0.43059185758737006,
        0.09805328339852282,
        0.8104023752957474,
        0.7673972319392628,
        0.1459688415802641,
        0.4529949031554191,
        0.4735878233525810,
        0.3312174135701585,
        0.8970506176961734,
    ];

    println!(
        "Sample point #1: {}",
        integrand.evaluate(&point1[..dims], 1.0, None,)
    );

    println!(
        "Sample point #2: {}",
        integrand.evaluate(&point2[..dims], 1.0, None,)
    );

    integrand.set_verbosity(0);

    let data = Process {
        integrand: integrand,
    };
    let r = ci.vegas(dims, 1, 1, CubaVerbosity::Progress, 0, data);

    println!("{:#?}", r);
}
