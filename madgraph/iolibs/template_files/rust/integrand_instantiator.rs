use std::collections::HashMap;
use std::env;

use crate::integrand::Integrand;
use crate::run_card::RunCard;
use crate::param_card::ParamCard;
use crate::settings_card::SettingsCard;

%(instantiate_integrands_header)s

pub struct AllIntegrands_%(output_name)s {
    pub all_integrands : HashMap<usize, Integrand>,
}

impl AllIntegrands_%(output_name)s {
    pub fn new() -> AllIntegrands_%(output_name)s {

        let settings_card = SettingsCard::new(&(env::var("PROC_ROOT").unwrap_or("..".to_string()).to_owned()+"/rust/Cards/settings.yaml"));
        let run_card = RunCard::new(&(settings_card.root_path.clone() + "/rust/Cards/run_card.yaml"));
        let param_card = ParamCard::new(&(settings_card.root_path.clone() + "/rust/Cards/param_card.yaml"));

        let param_card_path = settings_card.root_path.clone() + "/Cards/param_card.dat";

        let mut all_integrands = HashMap::new();

        %(instantiate_integrands)s

        AllIntegrands_%(output_name)s { all_integrands }

    }
}
