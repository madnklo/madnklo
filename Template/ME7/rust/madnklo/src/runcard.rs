use std::fs::File;
use std::io::BufReader;
use serde::Deserialize;
use serde_yaml;

#[derive(Deserialize)]
pub struct RunCard {
    draa: f64,
    draj: f64,
    dral: f64,
    drjj: f64,
    drjl: f64,
    drll: f64,
    dsqrt_q2fact1: f64,
    dsqrt_q2fact2: f64,
    dynamical_scale_choice: isize,
    ebeam1: f64,
    ebeam2: f64,
    etaa: f64,
    etaj: f64,
    etal: f64,
    fixed_fac_scale: bool,
    fixed_ren_scale: bool,
    flavor_cuts: Option<String>,
    fo_analysis: Option<String>,
    integrator: String,
    iseed: usize,
    lhaid: usize,
    lpp1: isize,
    lpp2: isize,
    maxjetflavor: usize,
    nevents: isize,
    nhel: usize,
    pdlabel: String,
    pta: f64,
    ptj: f64,
    ptl: f64,
    run_tag: String,
    scale: f64,
    use_syst: bool,
   
}

impl RunCard {
    pub fn new(filename: &str) -> RunCard {
        let f = File::open(filename).expect("Could not open run card");
        let reader = BufReader::new(f);
        serde_yaml::from_reader(reader).expect("Could not read run card")
    }
}