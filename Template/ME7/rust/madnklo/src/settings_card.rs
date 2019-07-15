use std::fs::File;
use std::io::BufReader;
use serde::Deserialize;
use serde_yaml;

#[derive(Deserialize)]
pub struct SettingsCard {
    pub root_path : String,
    pub lhapdf_library_path : String,
}

impl SettingsCard {
    pub fn new(filename: &str) -> SettingsCard {
        let f = File::open(filename).expect("Could not open the settings card");
        let reader = BufReader::new(f);
        serde_yaml::from_reader(reader).expect("Could not read the settings card")
    }
}