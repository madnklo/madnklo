use crate::param_card::ParamCard;
use crate::phase_space_generator::{FlatPhaseSpaceGenerator, PhaseSpaceGenerator};
use crate::run_card::RunCard;
use crate::settings_card::SettingsCard;
use epsilon_expansion::EpsilonExpansion;
use libc::{c_double, c_int};
use std::collections::HashMap;
use std::f64::consts::PI;
use std::ffi::CString;
use std::mem;
use typenum::{N6, U13};
use vector::LorentzVector;

mod LHAPDF {
    use libc::{c_char, c_double, c_int, c_void};

    #[link(name = "LHAPDF", kind = "static")]
    extern "C" {
        pub fn initpdfsetbyname_(setname: *const c_char, setnamelength: c_int);
        pub fn initpdf_(index: *const c_int);
        pub fn alphaspdf_(scale: *const c_double) -> f64;
        pub fn evolvepdf_(x: *const c_double, q: *const c_double, f: *mut c_double) -> f64;
        pub fn evolvepdfphoton_(
            x: *const c_double,
            q: *const c_double,
            f: *mut c_double,
            photon: *mut c_double,
        ) -> f64;
        pub fn has_photon_() -> bool;
    }
}

mod FJCORE {
    use libc::{c_char, c_double, c_int, c_void};

    #[link(name = "fjcore", kind = "static")]
    extern "C" {
        pub fn fjcoreppgenkt_(
            p: *const c_double,
            npart: *const c_int,
            R: *const c_double,
            palg: *const c_double,
            jets: *mut c_double,
            njets: *mut c_int,
        );
        pub fn fjcoreeegenkt_(
            p: *const c_double,
            npart: *const c_int,
            R: *const c_double,
            palg: *const c_double,
            jets: *mut c_double,
            njets: *mut c_int,
        );
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
    bjorken_x_rescalings: (Option<f64>, Option<f64>),
    is_a_mirrored_event: bool,
}

impl Event {
    /// Convolute the weights_per_flavor_configurations.
    fn apply_PDF_convolution(&mut self, mu_f: (f64, f64)) {
        // In order for the + distributions of the PDF counterterms and integrated
        // collinear ISR counterterms to hit the PDF only (and not the matrix elements or
        // observables functions), a change of variable is necessary: xb_1' = xb_1 * xi1
        // In this case, the xi1 to factor to apply is given by the Bjorken_x_rescalings
        // factor, which must be used both for rescaling the argument of the PDF as well
        // as bringing a factor 1/xi for the jacobian of the change of variable.
        for (flavors, flavour_weight) in self.weights_per_flavor_configurations.iter_mut() {
            let mut PDFs = 1.;

            PDFs *= Event::get_pdfQ(
                flavors.0[0],
                self.bjorken_xs.0 * self.bjorken_x_rescalings.0.unwrap(),
                mu_f.0,
            );
            PDFs *= self.bjorken_x_rescalings.0.unwrap();

            if flavors.0.len() == 2 {
                PDFs *= Event::get_pdfQ(
                    flavors.0[1],
                    self.bjorken_xs.1 * self.bjorken_x_rescalings.1.unwrap(),
                    mu_f.1,
                );
                PDFs *= self.bjorken_x_rescalings.1.unwrap();
            }

            *flavour_weight *= PDFs;
        }
    }

    /// Call the PDF and return the corresponding density.
    fn get_pdfQ(pdg: isize, x: f64, scale: f64) -> f64 {
        if pdg != 21 && pdg != 22 && (pdg == 0 || pdg.abs() > 6) {
            return 1.;
        }

        let mut f: [f64; 13] = [0.; 13]; // TODO: cache for all flavours
        unsafe {
            let pdf_weight = if LHAPDF::has_photon_() {
                let mut photon: f64 = 0.;
                LHAPDF::evolvepdfphoton_(
                    &x as *const c_double,
                    &scale as *const c_double,
                    &mut f[0] as *mut c_double,
                    &mut photon as *mut c_double,
                );
                match pdg {
                    22 => photon,
                    21 => f[6],
                    a => f[(a + 6) as usize],
                }
            } else {
                LHAPDF::evolvepdf_(
                    &x as *const c_double,
                    &scale as *const c_double,
                    &mut f[0] as *mut c_double,
                );
                match pdg {
                    22 => unreachable!("Selected pdf set does not include a photon"),
                    21 => f[6],
                    a => f[(a + 6) as usize],
                }
            };

            pdf_weight / x
        }
    }
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

    n_initial: usize,
    n_final: usize,
    n_unresolved_particles: usize,
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

    // level of debug info
    verbosity: usize,
}

impl Integrand {
    pub fn new(
        n_processes: usize,
        all_flavor_configurations: HashMap<usize, Vec<(Vec<isize>, Vec<isize>)>>,
        n_initial: usize,
        n_final: usize,
        n_unresolved_particles: usize, // TODO: introduce full contribution_definition
        masses: (Vec<f64>, Vec<f64>),
        run_card: RunCard,
        param_card: ParamCard,
        settings_card: SettingsCard,
        integrand_evaluator: Box<IntegrandEvaluator>,
    ) -> Integrand {
        // initialise the PDF, hardcoded for now
        unsafe {
            // note that if the file does not exist, lhapdf will crash without a clear error
            let name = CString::new("PDF4LHC15_nlo_30").unwrap();
            LHAPDF::initpdfsetbyname_(name.as_ptr(), name.to_bytes().len() as c_int);
            let v = 0 as c_int;
            LHAPDF::initpdf_(&v as *const c_int);
        }

        let collider_energy = run_card.ebeam1 + run_card.ebeam2;

        Integrand {
            n_processes,
            all_flavor_configurations,
            n_initial,
            n_final,
            n_unresolved_particles,
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
            verbosity: 0,
        }
    }

    pub fn set_verbosity(&mut self, verbosity: usize) {
        self.verbosity = verbosity;
    }

    pub fn get_dimensions(&self) -> usize {
        self.phase_space_generator.get_dimensions()
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

    ///    Implementation of a minimal set of isolation cuts. This can be made much nicer in the future and
    ///    will probably be taken outside of this class so as to ease user-specific cuts, fastjet, etc...
    ///    We consider here a two-level cuts system, this first one of which is flavour blind.
    ///    The 'n_jets_allowed_to_be_clustered' is an option that allows to overwrite the
    ///    maximum number of jets that can be clustered and which is by default taken to be:
    ///        self.contribution_definition.n_unresolved_particles
    ///    This is useful when using this function to apply cuts to the reduced PS of the CTs.
    ///    Finally, the options 'xb_1' and 'xb_2' allow to specify the boost bringing
    ///    the PS point back from the c.o.m. to the lab frame, necessary if some cuts are not
    ///    invariant under a boost along the beam axis.
    pub fn pass_flavor_blind_cuts(
        &self,
        PS_point: &[LorentzVector<f64>],
        process_pdg_final: &[isize],
        n_jets_allowed_to_be_clustered: usize,
        xb_1: Option<f64>,
        xb_2: Option<f64>,
    ) -> bool {
        // These cuts are not allowed to resolve flavour, but only whether a particle is a jet or not
        let is_a_jet = |pdg: isize| {
            let pdga = pdg.abs() as usize;
            pdga == 21 || (pdga >= 1 && pdga <= self.run_card.maxjetflavor)
        };

        fn is_a_lepton(pdg: isize) -> bool {
            let pdga = pdg.abs();
            pdga == 11 || pdga == 13 || pdga == 15
        }

        fn is_a_neutrino(pdg: isize) -> bool {
            let pdga = pdg.abs();
            pdga == 12 || pdga == 14 || pdga == 16
        }

        fn is_a_photon(pdg: isize) -> bool {
            pdg == 22
        }

        // If the cuts depend on the boost to the lab frame in case of hadronic collision
        // then the boost below can be used. Notice that the jet cuts is typically formulated in terms of
        // *pseudo* rapidities and not rapidities, so that when jets are clustered into a massive combined
        // momentum, this makes it important to boost to the lab frame. We therefore turn on this boost here
        // by default.
        if (xb_1 != None) && (xb_2 != None) && (xb_1, xb_2) != (Some(1.), Some(1.)) {
            panic!("Boosting to lab frame not supported yet");
            // Prevent border effects
            //PS_point = PS_point.get_copy()
            //PS_point.boost_from_com_to_lab_frame(xb_1, xb_2, self.run_card['ebeam1'], self.run_card['ebeam2'])
        }

        for (i, p) in PS_point[self.n_initial..].iter().enumerate() {
            if is_a_photon(process_pdg_final[i]) {
                if p.pt() < 100.0 {
                    return false;;
                }
            }
        }

        for (i, p) in PS_point[self.n_initial..].iter().enumerate() {
            if is_a_photon(process_pdg_final[i]) {
                for (j, p2) in PS_point[self.n_initial..].iter().enumerate() {
                    if j != i && process_pdg_final[j] != 21 {
                        if p.deltaR(p2) < 0.4 {
                            return false;
                        }
                    }
                }
            }
        }

        // JET CLUSTERING AND CUTS

        let ptj_cut = self.run_card.ptj;
        let drjj_cut = self.run_card.drjj;
        let etaj_cut = self.run_card.etaj;

        if ptj_cut <= 0. && drjj_cut <= 0. && etaj_cut <= 0. {
            return true;
        }

        let mut all_jets;
        if drjj_cut > 0. {
            // First identify all partonic jets
            let mut jets_list = vec![];
            for (i, p) in PS_point[self.n_initial..].iter().enumerate() {
                if is_a_jet(process_pdg_final[i]) {
                    jets_list.extend(&[p.t, p.x, p.y, p.z]);
                }
            }
            // Count partonic jets
            let starting_n_jets = jets_list.len();

            // Cluster them with fastjet
            let mut jets_flat = vec![0.; jets_list.len()];
            all_jets = Vec::with_capacity(jets_list.len() / 4);
            let mut actual_len: c_int = 0;
            let len: c_int = (jets_list.len() / 4) as c_int;
            let palg = -1.0;

            unsafe {
                FJCORE::fjcoreeegenkt_(
                    &jets_list[0] as *const c_double,
                    &len as *const c_int,
                    &drjj_cut as *const c_double,
                    &palg as *const c_double,
                    &mut jets_flat[0] as *mut c_double,
                    &mut actual_len as *mut c_int,
                );
            }

            let mut jet_count = 0;
            for i in 0..actual_len as usize {
                let jet = LorentzVector::from_slice(&jets_flat[i * 4..(i + 1) * 4]);
                // Remove jets whose pT is below the user defined cut
                if ptj_cut <= 0. || jet.pt() >= ptj_cut {
                    all_jets.push(jet);
                }
            }

            // Make sure that the number of clustered jets is at least larger or equal to the
            // starting list of jets minus the number of particles that are allowed to go
            // unresolved in this contribution.
            if all_jets.len() < starting_n_jets - n_jets_allowed_to_be_clustered {
                return false;
            }
        } else {
            all_jets = PS_point[self.n_initial..]
                .iter()
                .enumerate()
                .filter_map(|(i, p)| {
                    if is_a_jet(process_pdg_final[i]) {
                        Some(*p)
                    } else {
                        None
                    }
                })
                .collect();

            if ptj_cut > 0. {
                // Apply the Ptj cut first
                for p in &all_jets {
                    if p.pt() < ptj_cut {
                        return false;
                    }
                }
            }

            if drjj_cut > 0. {
                // And then the drjj cut
                for (i, p1) in all_jets.iter().enumerate() {
                    for p2 in &all_jets[i + 1..] {
                        if p1.deltaR(p2) < drjj_cut {
                            return false;
                        }
                    }
                }
            }
        }

        // Now handle all other cuts
        if etaj_cut > 0. {
            for p_jet in &all_jets {
                if p_jet.pseudoRap().abs() > etaj_cut {
                    return false;
                }
            }
        }

        // LEPTON AND PHOTON CUTS

        for (i, p) in PS_point[self.n_initial..].iter().enumerate() {
            // photons
            if is_a_photon(process_pdg_final[i]) {
                if self.run_card.pta > 0.0 && p.pt() < self.run_card.pta {
                    return false;
                }
                if self.run_card.etaa > 0.0 && p.pseudoRap().abs() > self.run_card.etaa {
                    return false;
                }
                for (j, p2) in PS_point[self.n_initial..].iter().enumerate() {
                    if j > i && is_a_photon(process_pdg_final[j]) {
                        if self.run_card.draa > 0.0 && p.deltaR(p2) < self.run_card.draa {
                            return false;
                        }
                    }
                    if is_a_lepton(process_pdg_final[j]) {
                        if self.run_card.dral > 0.0 && p.deltaR(p2) < self.run_card.dral {
                            return false;
                        }
                    }
                }
                for (j, p_jet) in all_jets.iter().enumerate() {
                    if self.run_card.draj > 0.0 && p.deltaR(p_jet) < self.run_card.draj {
                        return false;
                    }
                }
            }

            // leptons
            if is_a_lepton(process_pdg_final[i]) {
                if self.run_card.ptl > 0.0 && p.pt() < self.run_card.ptl {
                    return false;
                }
                if self.run_card.etal > 0.0 && p.pseudoRap().abs() > self.run_card.etal {
                    return false;
                }
                for (j, p2) in PS_point[self.n_initial..].iter().enumerate() {
                    if j <= i {
                        continue;
                    }
                    if is_a_lepton(process_pdg_final[j]) {
                        if self.run_card.drll > 0.0 && p.deltaR(p2) < self.run_card.drll {
                            return false;
                        }
                    }
                }
                for (j, p_jet) in all_jets.iter().enumerate() {
                    if self.run_card.drjl > 0.0 && p.deltaR(p_jet) < self.run_card.drjl {
                        return false;
                    }
                }
            }
        }

        // All cuts pass, therefore return true
        true
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
        x1: Option<f64>,
        x2: Option<f64>,
        sector_info: Option<Sector>,
    ) -> Vec<Event> {
        let mut sigma_wgt = base_weight;

        let alpha_s = self.param_card.aS;
        let all_flavor_configurations = &self.all_flavor_configurations[&process_id];

        // Apply flavor blind cuts
        if !self.pass_flavor_blind_cuts(
            &self.external_momenta,
            &all_flavor_configurations[0].1,
            self.n_unresolved_particles,
            None,
            None,
        ) {
            return vec![];
        }

        let ME_evaluation = self.integrand_evaluator.call_ME_for_key(
            &self.external_momenta_map,
            alpha_s,
            mu_r,
            process_id,
            0,
            0,
        );
        sigma_wgt *= ME_evaluation[0];

        if self.verbosity > 0 {
            println!("sigma={:e}", ME_evaluation[0]);
        }

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
            bjorken_x_rescalings: (x1, x2),
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

        let (PS_weight, x1s, x2s) = self
            .phase_space_generator
            .get_PS_point(random_variables, &mut self.external_momenta);

        for (i, x) in self.external_momenta.iter().enumerate() {
            self.external_momenta_map.insert(i + 1, *x);

            if self.verbosity > 0 {
                println!("PS_points {}={:e}", i + 1, x);
            }
        }

        wgt *= PS_weight;

        if self.verbosity > 0 {
            println!("PS_weight={:e}", PS_weight);
        }

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

        let alpha_s = unsafe { LHAPDF::alphaspdf_(&mu_r as *const c_double) };
        if self.verbosity > 0 {
            println!("alpha_s = {}", alpha_s);
        }

        // overwrite the parameters
        self.param_card.aS = alpha_s;
        //self.param_card.MU_R = mu_r; // FIXME: this key does not exist. Generate with Option<f64> instead!

        // Now loop over processes
        let mut total_wgt = 0.;
        for process_id in 0..self.n_processes {
            // Make sure to skip processes that should not be considered if the process_key is specified
            // FIXME: this part is not clear to me.
            let mut sector_info = None;
            /*if self.selected_sectors.len() == 1[process_id] let Some(selected_process_ids_for_sector) = process_ids_for_sector {
                if !process_ids_for_sector.contains(process_id) {
                    continue;
                }
                sector_info = self.process_ids_for_sector[process_id];
            }*/

            let mut events = self.sigma(
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

            // Apply flavor blind cuts
            let mut filtered_events = Vec::with_capacity(events.len());
            for e in events {
                // TODO: which index to take here for all_flavor_configurations?
                // the cut is flavour blind, so it doesn't matter?
                if self.pass_flavor_blind_cuts(
                    &self.external_momenta,
                    &self.all_flavor_configurations[&process_id][0].1,
                    self.n_unresolved_particles,
                    None,
                    None,
                ) {
                    filtered_events.push(e);
                }
            }
            events = filtered_events;

            for event in &mut events {
                // And epsilon-expansion
                // Select particular terms of the EpsilonExpansion terms stored as weights
                //event.select_epsilon_expansion_term('finite')

                // Apply PDF convolution
                if self.run_card.lpp1 == self.run_card.lpp2 && self.run_card.lpp1 == 1 {
                    event.apply_PDF_convolution((mu_f1, mu_f2));
                }

                // Make sure Bjorken-x rescalings don't matter anymore
                event.bjorken_x_rescalings = (None, None);

                // Apply flavor sensitive cuts
                // if not events.filter_flavor_configurations(self.pass_flavor_sensitive_cuts):
                //    if __debug__: logger.debug('Events failed the flavour-sensitive generation-level cuts.')
                //    if __debug__: logger.debug(misc.bcolors.GREEN + 'Returning a weight of 0. for this integrand evaluation.' + misc.bcolors.ENDC)
                //    continue

                // Aggregate the total weight to be returned to the integrator
                total_wgt += event
                    .weights_per_flavor_configurations
                    .values()
                    .sum::<f64>();
            }

            // Select a unique flavor configuration for each event
            // Not so important, can be skipped.
            //events.select_a_flavor_configuration()

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