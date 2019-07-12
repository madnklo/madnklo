use std::collections::HashMap;
use std::env;

%(instantiate_integrands_header)s

struct AllIntegrands_%(output_name)s {
    pub all_integrands : HashMap<usize,Integrand>,
}

impl AllIntegrands_%(output_name)s {
    fn new() -> AllIntegrands_%(output_name)s {

        settings_card = SettingsCard::new(&(env::var("PROC_ROOT").unwrap().to_owned()+"/Cards/settings.yaml"));
        run_card = RunCard(&(settings_card.root_path+"/Cards/run_card.yaml"));
        param_card = RunCard(&(settings_card.root_path+"/Cards/param_card.yaml"));

        let mut all_integrands = HashMap<usize,Integrand>();

        %(instantiate_integrands)s

        AllIntegrands_%(output_name)s { all_integrands };

    }
}
