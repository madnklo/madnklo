use std::collections::HashMap;

struct AllIntegrands_%(output_name)s {
    pub all_integrands : HashMap<usize,Integrand>,
}

impl AllIntegrands_%(output_name)s {
    fn new() -> AllIntegrands_%(output_name)s {
        let all_integrands = AllIntegrands_%(output_name)s { HashMap() };

        %(instantiate_integrands)s

        all_integrands
    }
}
