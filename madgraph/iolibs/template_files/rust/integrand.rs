
struct Integrand:

    // Process map. Simplified version since most of the corresponding information
    // extracted from the abstract process instance will be hard-coded here.
    int n_processes = %(n_processes)d,
    // Hard-coded information obtained for each process_id
    // TO BE HARD-CODED
    all_flavor_configurations = { process_id : [(initial_PDGs, final_PDGs), ...] }

    // List of sector IDs for each process (may change between runs)
    // It is important that this can be modified upon instantiation of the sectors.
    processes_per_sector = Option < [ { process_id : Sector }, ] >
    selected_sectors = None

    // One may want to chose different PS generators, so ideally also set upon
    // instantiation but it can be hardcoded for now
    n_initial = %(n_initial)d;
    n_final = %(n_final)d
    masses = []
    phase_space_generator = None

    // PDF object given by LHPADF. Will be None for e+ e-.
    pdf = None

    // Run card storing run information, Go through yaml maybe?
    run_card = RunCard()

impl Integrand {

    fn new(
        &self,
        processes_per_sector: Option < [ { process_id : Sector }, ] >,
        selected_sectors : Option< [sector_ids,] >,
        phase_space_generator: Option <PSGenerator>,
        ) {
            match processes_per_sector {
                Some(sectors_specifier) => {self.processes_per_sector = sectors_specifier};
                _ => {
                    self.processes_per_sector = [ { id : None for id in [1..(n_processes+1)]}];
                },
            };

            match phase_space_generator {
                Some(phase_space_generator) => {self.phase_space_generator = phase_space_generator};
                _ => {
                    // instantiate here a default flat PS generator
                    self.phase_space_generator = FlatPSGenerate(...),
                },
            };

            self.selected_sectors = selected_sectors

            // Instantiate a run_card struct being the rust equivalent of run_card.dat
            self.run_card = ...

            // Get here a PDF instance from LHAPDF (use run_card to know which one to use).
            self.pdf = ...

        }

    impl call(&self, random_variables: f64[ %(n_random_variables)d], integrator_weight: Option < f64 > ) -> f64:

        integrator_jacobian = 1.
        match integrator_weight {
            Some(wgt) => {integrator_jacobian = wgt};
            _ => {},
        };

        final_weight = 0.
        if let Some(selected_sectors):
            final_weight += self.evaluate(random_variables, integrator_jacobian, None)
        else:
            for i_sector in self.selected_sectors:
                final_weight += self.evaluate(random_variables, integrator_jacobian, processes_per_sector[i_sector])


    impl evaluate(
        random_variables: f64[ %(n_random_variables)d],
        integrator_weight: f64,
        selected_process_ids_for_sector: Option < { process_id : Sector } >,
    ) -> f64:

        // A unique float must be returned
        wgt = 1.0
        // And the conversion from GeV^-2 to picobarns
        wgt *= 0.389379304e9

        // Below needs to expand the applicability of the PS generator and export the beam
        // information as well.
        // WARNING: We must generate here a PS_point which is a dictionary with integer IDs
        PS_point, PS_weight, x1s, x2s = self.phase_space_generator.get_PS_point(PS_random_variables)

        // Flux: trivial taken from masses information above

        // Need to implement a token get_scales function taking inputs from self.run_card
        mu_r, mu_f1, mu_f2 = self.get_scales(PS_point)

        // Code alpha_s_runner function using LHAPDF if enabled (see later on how to link LHAPDF)
        alpha_s = self.alpha_s_runner(mu_r)

        // Now loop over processes
        total_wgt = 0.
        for process_id in [1..(self.n_processes+1)]:

            // Make sure to skip processes that should not be considered if the process_key is specified
            sector_info = None
            if let process_ids_for_sector = Some(selected_process_ids_for_sector):
                if process_id not in process_ids_for_sector:
                    continue
                sector_info = process_ids_for_sector[process_id]

            all_flavor_configurations = self.all_flavor_configurations[process_id]

            // This sigma function for Born will mostly have to do with finding a way to
            // hardcode the call to the ME
            // IMPORTANT to build an Event rust equivalent
            events = self.sigma(
                PS_point, process_id all_flavor_configurations,
                wgt, mu_r, mu_f1, mu_f2, xb_1, xb_2, xi1, xi2, sector_info)

            if len(events)==0:
                continue
            // Now handle the process mirroring, by adding, for each ME7Event defined with
            // requires_mirroring = True, a new mirrored event with the initial states swapped,
            // as well as the Bjorken x's and rescalings and with p_z -> -p_z on all final
            // state momenta
            events.generate_mirrored_events()

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
            events.select_epsilon_expansion_term('finite')

            // Apply PDF convolution
            // For this one must use the LHPADF library whose location is given by MG5aMC options
            // and for this you must also use the "run_card" Rust run_card object
            if self.run_card['lpp1']==self.run_card['lpp2']==1:
                events.apply_PDF_convolution( self.get_pdfQ2,
                                               (self.pdf, self.pdf), (mu_f1**2, mu_f2**2) )

            // Make sure Bjorken-x rescalings don't matter anymore
            // Simple safety feature to add
            for event in events:
                event.set_Bjorken_rescalings(None, None)

            // Apply flavor sensitive cuts
            // if not events.filter_flavor_configurations(self.pass_flavor_sensitive_cuts):
            //    if __debug__: logger.debug('Events failed the flavour-sensitive generation-level cuts.')
            //    if __debug__: logger.debug(misc.bcolors.GREEN + 'Returning a weight of 0. for this integrand evaluation.' + misc.bcolors.ENDC)
            //    continue

            // Aggregate the total weight to be returned to the integrator
            total_wgt += events.get_total_weight()

            // Select a unique flavor configuration for each event
            // Not so important, can be skipped.
            //events.select_a_flavor_configuration()

            // For debugging it would be useful to have the possibility of monitoring all events
            // generated
            if __debug__: logger.debug('Short-distance events for this subprocess:\n%s\n'%str(events))

            // Finally apply observables in a thread-safe manner. The only really safe way will
            // probabaly be to write to different files that are later recombined and not use locks
            // if self.apply_observables:
            //    events.apply_observables(self.observable_list, integrator_jacobian)

        // Useful to be able to debug final weight returned for this contribution
        // if __debug__: logger.debug(misc.bcolors.GREEN + "Final weight returned: %.5e"%total_wgt + misc.bcolors.ENDC)

        return total_wgt
}


