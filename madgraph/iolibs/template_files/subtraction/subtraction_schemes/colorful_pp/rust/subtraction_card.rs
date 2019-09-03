use std::fs::File;
use std::io::BufReader;
use serde::Deserialize;
use serde_yaml;

#[derive(Deserialize)]
pub struct SubtractionCard {
    pub alpha_0: f64,
    pub y_0: f64,
    pub y_0_prime: f64,
    pub d_0: f64,
    pub d_0_prime: f64,
}

impl SubtractionCard {
    pub fn new(filename: &str) -> SubtractionCard {
        let f = File::open(filename).expect("Could not open the subtraction card");
        let reader = BufReader::new(f);
        serde_yaml::from_reader(reader).expect("Could not read the subtraction card")
    }
}