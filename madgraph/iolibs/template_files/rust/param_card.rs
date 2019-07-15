use std::fs::File;
use std::io::BufReader;
use serde::Deserialize;
use serde_yaml;

#[derive(Deserialize)]
pub struct ParamCard {
%(model_parameters)s
}

impl ParamCard {
    pub fn new(filename: &str) -> ParamCard {
        let f = File::open(filename).expect("Could not open run card");
        let reader = BufReader::new(f);
        serde_yaml::from_reader(reader).expect("Could not read run card")
    }
}