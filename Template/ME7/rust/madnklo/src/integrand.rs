use crate::param_card::ParamCard;
use crate::phase_space_generator::{FlatPhaseSpaceGenerator, PhaseSpaceGenerator};
use crate::run_card::RunCard;
use crate::settings_card::SettingsCard;
use epsilon_expansion::EpsilonExpansion;
use std::collections::HashMap;
use std::f64::consts::PI;
use std::ffi::CString;
use std::mem;
use typenum::{N6, U13};
use vector::LorentzVector;

mod LHAPDF {
    use libc::{c_char, c_double, c_int, c_longlong, c_void};

    #[link(name = "LHAPDF")]
    extern "C" {
        pub fn initpdfsetbyname(p: *const c_char);
        pub fn initpdf(index: c_int);
        pub fn alphasQ(scale: f64) -> f64;
    }
}

pub trait IntegrandEvaluator {
    fn call_ME_for_key(
        &mut self,
        p: &HashMap<usize, LorentzVector<f64>>,
        alpha_s: f64,
        mu_r: f64,
        i_process: usize,
        i_CT: usize,
        i_call: usize,
    ) -> EpsilonExpansion<U13, N6>;
}

pub struct Sector {
    identifier: String,
    n_legs: usize,
}

pub struct Event {
    ps_point: HashMap<usize, LorentzVector<f64>>,
    weights_per_flavor_configurations: HashMap<(Vec<isize>, Vec<isize>), f64>, // FIXME: a terrible type!
    host_contribution_definition: String,
    requires_mirroring: bool,
    bjorken_xs: (f64, f64),
    bjorken_x_rescalings: (f64, f64),
    is_a_mirrored_event: bool,
}

pub struct Integrand {
    // Process map. Simplified version since most of the corresponding information
    // extracted from the abstract process instance will be hard-coded here.
    n_processes: usize,
    // Information obtained for each process_id
    all_flavor_configurations: HashMap<usize, Vec<(Vec<isize>, Vec<isize>)>>,

    // List of sector IDs for each process (may change between runs)
    // It is important that this can be modified upon instantiation of the sectors.
    processes_per_sector: Vec<HashMap<usize, Option<Sector>>>,
    selected_sectors: Vec<usize>,

    // One may want to chose different PS generators, so ideally also set upon
    // instantiation but it can be hardcoded for now
    n_initial: usize,
    n_final: usize,
    masses: (Vec<f64>, Vec<f64>),
    phase_space_generator: Box<PhaseSpaceGenerator>,
    external_momenta: Vec<LorentzVector<f64>>,
    external_momenta_map: HashMap<usize, LorentzVector<f64>>,

    collider_energy: f64,

    // Cards storing run, model and settings information
    run_card: RunCard,
    param_card: ParamCard,
    settings_card: SettingsCard,

    // Integrand evaluator containing all the hardcoded information
    integrand_evaluator: Box<IntegrandEvaluator>,
}

impl Integrand {
    pub fn new(
        n_processes: usize,
        all_flavor_configurations: HashMap<usize, Vec<(Vec<isize>, Vec<isize>)>>,
        n_initial: usize,
        n_final: usize,
        masses: (Vec<f64>, Vec<f64>),
        run_card: RunCard,
        param_card: ParamCard,
        settings_card: SettingsCard,
        integrand_evaluator: Box<IntegrandEvaluator>,
    ) -> Integrand {
        // initialise the PDF, hardcoded for now
        unsafe {
            let name = CString::new("PDF4LHC15_nlo_30").unwrap();
            LHAPDF::initpdfsetbyname(name.as_ptr());
            LHAPDF::initpdf(0);
        }

        let collider_energy = run_card.ebeam1 + run_card.ebeam2;

        Integrand {
            n_processes,
            all_flavor_configurations,
            n_initial,
            n_final,
            masses: masses.clone(),
            run_card,
            param_card,
            settings_card,
            collider_energy,
            // TODO: take values from run card or param card?
            phase_space_generator: Box::new(FlatPhaseSpaceGenerator::new(
                n_initial,
                masses,
                collider_energy,
                (0, 0),
                false,
                (false, false),
            )),
            external_momenta: vec![LorentzVector::default(); n_initial + n_final],
            external_momenta_map: (1..=n_initial + n_final)
                .map(|i| (i, LorentzVector::default()))
                .collect(),
            processes_per_sector: vec![(0..n_processes).map(|id| (id, None)).collect()],
            selected_sectors: vec![],
            integrand_evaluator,
        }
    }

    pub fn set_sectors(
        &mut self,
        processes_per_sector: Vec<HashMap<usize, Option<Sector>>>,
        selected_sectors: Vec<usize>,
    ) {
        self.processes_per_sector = processes_per_sector;
        self.selected_sectors = selected_sectors;
    }

    pub fn set_ps_generator(&mut self, phase_space_generator: Box<PhaseSpaceGenerator>) {
        self.phase_space_generator = phase_space_generator;
    }

    pub fn call(&mut self, random_variables: &[f64], integrator_weight: Option<f64>) -> f64 {
        let integrator_jacobian = match integrator_weight {
            Some(wgt) => wgt,
            _ => 1.,
        };

        let mut final_weight = 0.;
        if self.selected_sectors.is_empty() {
            final_weight = self.evaluate(random_variables, integrator_jacobian, None);
        } else {
            let mut sectors = mem::replace(&mut self.selected_sectors, vec![]);
            for i_sector in &sectors {
                final_weight +=
                    self.evaluate(random_variables, integrator_jacobian, Some(*i_sector));
            }
            mem::swap(&mut self.selected_sectors, &mut sectors);
        }

        final_weight
    }

    #[inline]
    /// Kahlen function.
    fn lambda(s: f64, sqrMA: f64, sqrMB: f64) -> f64 {
        s * s + sqrMA * sqrMA + sqrMB * sqrMB - 2. * s * sqrMA - 2. * sqrMB * sqrMA - 2. * s * sqrMB
    }

    pub fn sigma(
        &mut self,
        process_id: usize,
        base_weight: f64,
        mu_r: f64,
        mu_f1: f64,
        mu_f2: f64,
        xb_1: f64,
        xb_2: f64,
        x1: f64,
        x2: f64,
        sector_info: Option<Sector>,
    ) -> Vec<Event> {
        let mut sigma_wgt = base_weight;

        let alpha_s = self.param_card.aS;
        let all_flavor_configurations = &self.all_flavor_configurations[&process_id];

        // Apply flavor blind cuts
        //if not self.pass_flavor_blind_cuts(PS_point, all_flavor_configurations[0]):
        //    if __debug__: logger.debug('Event failed the flavour_blind generation-level cuts.')
        //    if __debug__: logger.debug(misc.bcolors.GREEN + 'Returning a weight of 0. for this integrand evaluation.' + misc.bcolors.ENDC)
        //    return None

        let ME_evaluation = self.integrand_evaluator.call_ME_for_key(
            &self.external_momenta_map,
            alpha_s,
            mu_r,
            process_id,
            0,
            0,
        );
        sigma_wgt *= ME_evaluation[0];

        // Return the lone LO event
        vec![Event {
            ps_point: self.external_momenta_map.clone(),
            host_contribution_definition: "N/A".to_owned(),
            weights_per_flavor_configurations: all_flavor_configurations
                .iter()
                .map(|f| (f.clone(), sigma_wgt))
                .collect(),
            requires_mirroring: false, // FIXME: hardcoding to false
            bjorken_xs: (xb_1, xb_2),
            bjorken_x_rescalings: (xb_1, xb_2), // FIXME: what to do here?
            is_a_mirrored_event: false,
        }]
    }

    pub fn evaluate(
        &mut self,
        random_variables: &[f64],
        integrator_weight: f64,
        selected_process_ids_for_sector: Option<usize>,
    ) -> f64 {
        // Return a weight in picobarns (from GeV^-2)
        let mut wgt = 0.389379304e9;

        // Below needs to expand the applicability of the PS generator and export the beam
        // information as well.
        // WARNING: We must generate here a PS_point which is a dictionary with integer IDs
        // TODO: add the integer ids later
        let (PS_weight, x1s, x2s) = self
            .phase_space_generator
            .get_PS_point(random_variables, &mut self.external_momenta);

        for (i, x) in self.external_momenta.iter().enumerate() {
            self.external_momenta_map.insert(i + 1, *x);
        }

        wgt *= PS_weight;

        let (xb_1, xi1) = x1s;
        let (xb_2, xi2) = x2s;

        // The E_cm entering the flux factor is computed *without* including the xi<i> rescalings
        // This is necessary and must be kept so.
        let E_cm = (xb_1 * xb_2).sqrt() * self.collider_energy;

        // Include the flux factor
        let mut flux = if self.n_initial == 2 {
            1. / (2.
                * (Integrand::lambda(
                    E_cm * E_cm,
                    self.masses.0[0] * self.masses.0[0],
                    self.masses.0[1] * self.masses.0[1],
                ))
                .sqrt())
        } else {
            1. / (2. * E_cm)
        };
        flux /= (2. * PI).powi(3 * self.n_final as i32 - 4);
        wgt *= flux;

        // get the fixed scale functions
        let (mu_r, mu_f1, mu_f2) = (
            self.run_card.scale,
            self.run_card.dsqrt_q2fact1,
            self.run_card.dsqrt_q2fact2,
        );

        let alpha_s = unsafe { LHAPDF::alphasQ(mu_r) };

        // overwrite the parameters
        self.param_card.aS = alpha_s;
        //self.param_card.MU_R = mu_r; // FIXME: this key does not exist. Generate with Option<f64> instead!

        // Now loop over processes
        let mut total_wgt = 0.;
        for process_id in 1..=self.n_processes {
            // Make sure to skip processes that should not be considered if the process_key is specified
            // FIXME: this part is not clear to me.
            let mut sector_info = None;
            /*if self.selected_sectors.len() == 1[process_id] let Some(selected_process_ids_for_sector) = process_ids_for_sector {
                if !process_ids_for_sector.contains(process_id) {
                    continue;
                }
                sector_info = self.process_ids_for_sector[process_id];
            }*/

            // This sigma function for Born will mostly have to do with finding a way to
            // hardcode the call to the ME
            // IMPORTANT to build an Event rust equivalent
            let events = self.sigma(
                process_id,
                wgt,
                mu_r,
                mu_f1,
                mu_f2,
                xb_1,
                xb_2,
                xi1,
                xi2,
                sector_info,
            );

            if events.is_empty() {
                continue;
            }
            // Now handle the process mirroring, by adding, for each ME7Event defined with
            // requires_mirroring = True, a new mirrored event with the initial states swapped,
            // as well as the Bjorken x's and rescalings and with p_z -> -p_z on all final
            // state momenta
            //events.generate_mirrored_events()

            // For things that we don't want to immediatly handle, we should just put placeholders
            // Apply flavor blind cuts
            // TODO implement flavour cuts
            /*
            if not events.filter_with_flavor_blind_cuts(self.pass_flavor_blind_cuts, process_pdgs,
                n_jets_allowed_to_be_clustered  = self.contribution_definition.n_unresolved_particles):
                if __debug__: logger.debug('All events failed the flavour_blind generation-level cuts.')
                if __debug__: logger.debug(misc.bcolors.GREEN + 'Returning a weight of 0. for this integrand evaluation.' + misc.bcolors.ENDC)
                continue
            */

            // And epsilon-expansion
            // Select particular terms of the EpsilonExpansion terms stored as weights
            //events.select_epsilon_expansion_term('finite')

            // Apply PDF convolution
            // For this one must use the LHPADF library whose location is given by MG5aMC options
            // and for this you must also use the "run_card" Rust run_card object
            //if self.run_card.lpp1 == self.run_card.lpp2 && self.run_card.lpp1 == 1:
            //    events.apply_PDF_convolution( self.get_pdfQ2, (self.pdf, self.pdf), (mu_f1**2, mu_f2**2) )

            // Make sure Bjorken-x rescalings don't matter anymore
            // Simple safety feature to add
            //for event in &mut events {
            //    event.set_Bjorken_rescalings(None, None)
            //}

            // Apply flavor sensitive cuts
            // if not events.filter_flavor_configurations(self.pass_flavor_sensitive_cuts):
            //    if __debug__: logger.debug('Events failed the flavour-sensitive generation-level cuts.')
            //    if __debug__: logger.debug(misc.bcolors.GREEN + 'Returning a weight of 0. for this integrand evaluation.' + misc.bcolors.ENDC)
            //    continue

            // Aggregate the total weight to be returned to the integrator
            total_wgt += events
                .iter()
                .map(|e| e.weights_per_flavor_configurations.values().sum::<f64>())
                .sum::<f64>();

            // Select a unique flavor configuration for each event
            // Not so important, can be skipped.
            //events.select_a_flavor_configuration()

            // For debugging it would be useful to have the possibility of monitoring all events
            // generated
            //if __debug__: logger.debug('Short-distance events for this subprocess:\n%s\n'%str(events))

            // Finally apply observables in a thread-safe manner. The only really safe way will
            // probabaly be to write to different files that are later recombined and not use locks
            // if self.apply_observables:
            //    events.apply_observables(self.observable_list, integrator_jacobian)
        }

        // Useful to be able to debug final weight returned for this contribution
        // if __debug__: logger.debug(misc.bcolors.GREEN + "Final weight returned: %.5e"%total_wgt + misc.bcolors.ENDC)
        total_wgt
    }
}
