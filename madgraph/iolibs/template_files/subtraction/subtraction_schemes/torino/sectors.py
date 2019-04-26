#
# Generic classes for the sector implementation in subtraction schemes
#

import commons.generic_sectors as generic_sectors

class Sector(generic_sectors.GenericSector):
    """ Class implementing a particular sector, with attributes identifying it and methods 
    for its evaluation."""

    def __init__(self, external_leg_numbers, **opts):
        super(Sector, self).__init__(**opts)
        self.external_leg_numbers = external_leg_numbers

    def __call__(self, PS_point, PDGs, counterterm_index=-1, input_mapping_index=-1):
        """ Given a kinematic configuration (PS_point is a *dict* here) and flavors of external states, returns the sector weight (float)
        from this sector function.
        Notice that the sectoring function will be called for both the real(s) event and *also* with the reduced
        kinematics of the counterevent, in which case it may be necessary to know which counterterm this reduced
        kinematics comes from, so that the counterterm_index is specified (not the full counterterm, because
        this function should eventually be coded at low-level and do minimal logic, so please do all the logic 
        the sector generator function. -1 for the counterterm_index implies no counterterm.
        The flavor mapping is set to None if a local counterterm is considered, and to the corresponding input
        mapping index if an integrated counterterm is considered.
        """

        sector_weight = 1.0

        return sector_weight

class SectorGenerator(generic_sectors.GenericSectorGenerator):
    """ Class responsible for generating the correct list of sectors to consider for specific processes."""

    def __init__(self, *args, **opts):
        super(SectorGenerator, self).__init__(*args, **opts)


    def __call__(self, contrib_definition, defining_process, counterterms=None, integrated_counterterms=None):
        """ Given a particular contribution definition, a particular defining process instance and the counterterms contained, this
        function must build the list of sectors to consider (or None if the subtraction
        does not require any) as well as the list of counterterms to consider for each of them.
        The list returned has the following format:

            all_sectors = [ sector_1, sector_2, ... ]

        where each sector_<i> is a dictionary with the following format:

            sector_<i> = {
                'sector'    :   sector_instance_or_identifier (exact format is up to the subtraction_scheme but must be callable)
                'counterterms' : [list of counterterm indices (as per the ordering in self.counterterms of the integrand) to be considered],
                'integrated_counterterms' : [list of tuples (counterterm indices, input_mapping_index) for integrated CT to be considered]
            }

        """
        
        model = defining_process.get('model')
        initial_state_PDGs, final_state_PDGs = defining_process.get_initial_final_ids()

        # Hardcoding to FKS-like sectors for e+(1) e-(2) > d(3) d~(4) g(5)
        all_sectors = []

        for sector_legs in [(5,3), (5,4)]:

            a_sector = {
                'sector' : None,
                'counterterms' : None,
                'integrated_counterterms' : None
            }

            a_sector['sector'] = Sector(external_leg_numbers=sector_legs)
            if counterterms is not None:
                a_sector['counterterms'] = range(len(counterterms))
            if integrated_counterterms is not None:
                a_sector['integrated_counterterms'] = [ (i_ct, i_mapping) for i_ct, ct in enumerate(integrated_counterterms)
                                                        for i_mapping in range(len(ct['input_mappings']))]

            all_sectors.append(a_sector)

        return all_sectors
