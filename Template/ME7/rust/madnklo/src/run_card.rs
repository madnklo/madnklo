use std::fs::File;
use std::io::BufReader;
use serde::Deserialize;
use serde_yaml;

#[derive(Deserialize)]
pub struct RunCard {
    pub draa: f64,
    pub draj: f64,
    pub dral: f64,
    pub drjj: f64,
    pub drjl: f64,
    pub drll: f64,
    pub dsqrt_q2fact1: f64,
    pub dsqrt_q2fact2: f64,
    pub dynamical_scale_choice: isize,
    pub ebeam1: f64,
    pub ebeam2: f64,
    pub etaa: f64,
    pub etaj: f64,
    pub etal: f64,
    pub fixed_fac_scale: bool,
    pub fixed_ren_scale: bool,
    pub flavor_cuts: Option<String>,
    pub fo_analysis: Option<String>,
    pub integrator: String,
    pub iseed: usize,
    pub lhaid: usize,
    pub lpp1: isize,
    pub lpp2: isize,
    pub maxjetflavor: usize,
    pub nevents: isize,
    pub nhel: usize,
    pub pdlabel: String,
    pub pta: f64,
    pub ptj: f64,
    pub ptl: f64,
    pub run_tag: String,
    pub scale: f64,
    pub use_syst: bool,
}

impl RunCard {
    pub fn new(filename: &str) -> RunCard {
        let f = File::open(filename).expect("Could not open run card");
        let reader = BufReader::new(f);
        serde_yaml::from_reader(reader).expect("Could not read run card")
    }
}