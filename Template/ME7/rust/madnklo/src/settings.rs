use std::fs::File;
use std::io::BufReader;
use serde::Deserialize;
use serde_yaml;

#[derive(Deserialize)]
pub struct SettingsCard {
    root_path : String,
    lhapdf_library_path : String,
}

impl SettingsCard {
    pub fn new(filename: &str) -> RunCard {
        let f = File::open(filename).expect("Could not open the settings card {}".format(filename));
        let reader = BufReader::new(f);
        serde_yaml::from_reader(reader).expect("Could not read the settings card {}".format(filename))
    }
}