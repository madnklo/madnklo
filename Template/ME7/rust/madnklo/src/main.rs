extern crate clap;
extern crate libc;
extern crate madnklo;
extern crate rand;

use cuba::{CubaIntegrator, CubaVerbosity};

use madnklo::all_integrands;
use madnklo::integrand::Integrand;

use clap::{App, Arg, ArgMatches, SubCommand};
use rand::prelude::*;
use std::str::FromStr;
use std::time::Instant;

#[derive(Debug, Clone)]
pub struct IntegratorSettings {
    pub integrator: Integrator,
    pub n_vec: usize,
    pub n_increase: usize,
    pub n_max: usize,
    pub n_start: usize,
    pub eps_rel: f64,
    pub eps_abs: f64,
    pub border: f64,
    pub maxpass: usize,
    pub maxchisq: f64,
    pub mindeviation: f64,
    pub n_new: usize,
    pub n_min: usize,
    pub flatness: f64,
    pub seed: i32,
    pub integrated_phase: IntegratedPhase,
    pub state_filename_prefix: Option<String>,
    pub survey_n_points: usize,
    pub survey_n_iterations: usize,
    pub refine_n_runs: usize,
    pub refine_n_points: usize,
    pub keep_state_file: bool,
    pub reset_vegas_integrator: bool,
    pub use_only_last_sample: bool,
}

#[derive(Debug, Clone)]
pub enum Integrator {
    Vegas,
    Suave,
    Cuhre,
    Divonne,
}

#[derive(Debug, Copy, Clone, PartialEq)]
pub enum IntegratedPhase {
    Real,
    Imag,
    Both,
}

impl Default for IntegratorSettings {
    fn default() -> IntegratorSettings {
        IntegratorSettings {
            integrator: Integrator::Vegas,
            n_increase: 0,
            n_vec: 1,
            n_start: 10000,
            n_max: 10000000,
            eps_rel: 1e-4,
            eps_abs: 0.,
            border: 1e-12,
            maxpass: 5,
            maxchisq: 10.,
            mindeviation: 0.25,
            n_new: 1000,
            n_min: 2,
            flatness: 50.,
            seed: 1,
            integrated_phase: IntegratedPhase::Real,
            state_filename_prefix: None,
            survey_n_points: 0,
            survey_n_iterations: 0,
            refine_n_runs: 0,
            refine_n_points: 0,
            keep_state_file: false,
            reset_vegas_integrator: true,
            use_only_last_sample: false,
        }
    }
}

struct Process<'a> {
    integrand: &'a mut Integrand,
}

#[inline(always)]
fn integrand_fn(
    x: &[f64],
    f: &mut [f64],
    process: &mut Process,
    _nvec: usize,
    _core: i32,
) -> Result<(), &'static str> {
    let res = process.integrand.evaluate(x, 1.0, None);

    if !res.is_finite() {
        println!("nan at x={:?}", x);
    }

    f[0] = res;

    Ok(())
}

/// Perform a single core benchmark.
fn bench(integrand: &mut Integrand, settings: &IntegratorSettings) {
    let mut x = vec![0.; integrand.get_dimensions()];
    let mut rng = rand::thread_rng();

    integrand.set_verbosity(0);

    let now = Instant::now();
    for _ in 0..settings.n_max {
        for xi in x.iter_mut() {
            *xi = rng.gen();
        }

        let _r = integrand.evaluate(&x, 1.0, None);
    }

    println!("{:#?}", now.elapsed());
}

/// Inspect a single point.
fn inspect<'a>(integrand: &mut Integrand, _settings: &IntegratorSettings, matches: &ArgMatches<'a>) {
    let x: Vec<_> = matches
        .values_of("point")
        .unwrap()
        .map(|x| f64::from_str(x).unwrap())
        .collect();
    if x.len() != integrand.get_dimensions() {
        panic!(
            "Dimension of the input point is incorrect. It should be {} but is {}.",
            integrand.get_dimensions(),
            x.len()
        );
    }

    integrand.set_verbosity(1);
    let r = integrand.evaluate(&x, 1.0, None);
    println!("Result={}", r);
}

fn main() {
    let matches = App::new("MadNkLO Rust backend")
        .version("0.1")
        .about("Numerically integrate your favourite processes")
        .arg(
            Arg::with_name("cores")
                .short("c")
                .long("cores")
                .value_name("NUMCORES")
                .help("Set the number of cores"),
        )
        .arg(
            Arg::with_name("samples")
                .short("s")
                .long("samples")
                .value_name("SAMPLES")
                .help("Number of samples per iteration"),
        )
        .arg(
            Arg::with_name("nstart")
                .long("nstart")
                .value_name("NSTART")
                .help("Number of starting samples per iteration"),
        )
        .arg(
            Arg::with_name("nincrease")
                .long("nincrease")
                .value_name("NINCREASE")
                .help("The increase in number of samples per iteration"),
        )
        .arg(
            Arg::with_name("seed")
                .long("seed")
                .value_name("SEED")
                .help("Specify the integration seed"),
        )
        .subcommand(SubCommand::with_name("bench").about("Run a benchmark"))
        .subcommand(
            SubCommand::with_name("inspect")
                .about("Inspect a single input point")
                .arg(Arg::with_name("point").required(true).min_values(3)),
        )
        .get_matches();

    let mut settings = IntegratorSettings::default();

    let mut cores = 0;
    if let Some(x) = matches.value_of("cores") {
        cores = usize::from_str(x).unwrap();
    }

    if let Some(x) = matches.value_of("seed") {
        settings.seed = i32::from_str(x).unwrap();
    }

    if let Some(x) = matches.value_of("samples") {
        settings.n_max = usize::from_str(x).unwrap();
    }

    if let Some(x) = matches.value_of("nstart") {
        settings.n_start = usize::from_str(x).unwrap();
    }

    if let Some(x) = matches.value_of("nincrease") {
        settings.n_increase = usize::from_str(x).unwrap();
    }

    let mut integrands = all_integrands::AllIntegrands_RUSTTEST::new();
    let integrand = integrands.all_integrands.get_mut(&1).unwrap();
    let dims = integrand.get_dimensions();

    if let Some(_) = matches.subcommand_matches("bench") {
        bench(integrand, &settings);
        return;
    }

    if let Some(matches) = matches.subcommand_matches("inspect") {
        inspect(integrand, &settings, matches);
        return;
    }

    let mut ci = CubaIntegrator::new(integrand_fn);
    ci.set_mineval(10)
        .set_nstart(settings.n_start as i64)
        .set_nincrease(settings.n_increase as i64)
        .set_maxeval(settings.n_max as i64)
        .set_epsrel(settings.eps_rel)
        .set_epsabs(settings.eps_abs)
        .set_border(settings.border)
        .set_maxpass(settings.maxpass as i32)
        .set_maxchisq(settings.maxchisq)
        .set_mindeviation(settings.mindeviation)
        .set_seed(settings.seed)
        .set_cores(cores, 1000);

    integrand.set_verbosity(0);

    let data = Process {
        integrand: integrand,
    };
    let r = ci.vegas(
        dims,
        if settings.integrated_phase == IntegratedPhase::Both {
            2
        } else {
            1
        },
        settings.n_vec,
        CubaVerbosity::Progress,
        0,
        data,
    );

    println!("{:#?}", r);
}
