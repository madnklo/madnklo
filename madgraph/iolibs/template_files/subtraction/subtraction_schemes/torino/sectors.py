#
# Generic classes for the sector implementation in subtraction schemes
#
import copy


import commons.generic_sectors as generic_sectors
import madgraph.various.misc as misc
import madgraph.core.diagram_generation as diagram_generation
import madgraph.fks.fks_common as fks_common
import madgraph.integrator.vectors as vectors
import logging


logger = logging.getLogger('madgraph')


def get_sector_wgt(q, p_sector):
    """ the sector function sigma_ij of the Torino paper (eq. 2.13-2.14)
     without the normalisation.
     - q is the total momentum of the incoming particles
     - p_sector is a list of the momenta of the particles defining the sector
    """
    s = q.square()
    s_ij = 2 * p_sector[0].dot(p_sector[1])
    s_qi = 2 * p_sector[0].dot(q)
    s_qj = 2 * p_sector[1].dot(q)
    e_i = s_qi / s
    w_ij = s * s_ij / s_qi / s_qj

    return 1 / e_i / w_ij


def get_sector_wgt_S(q, p_sector):
    """ the soft sector function sigma_ij of the Torino paper (eq. 2.15left)
     without the normalisation.
     - q is the total momentum of the incoming particles
     - p_sector is a list of the momenta of the particles defining the sector
    """
    s = q.square()
    s_ij = 2 * p_sector[0].dot(p_sector[1])
    s_qi = 2 * p_sector[0].dot(q)
    s_qj = 2 * p_sector[1].dot(q)
    w_ij = s * s_ij / s_qi / s_qj

    return 1 / w_ij


def get_sector_wgt_C(q, p_sector):
    """ the collinear sector function sigma_ij of the Torino paper (eq. 2.15right)
     without the normalisation.
     - q is the total momentum of the incoming particles
     - p_sector is a list of the momenta of the particles defining the sector
    """
    s = q.square()
    s_ij = 2 * p_sector[0].dot(p_sector[1])
    s_qi = 2 * p_sector[0].dot(q)
    s_qj = 2 * p_sector[1].dot(q)
    e_i = s_qi / s
    e_j = s_qj / s

    return 1. # instead of implementing 2.15 right, simply return 1 when summing i,j + j,i
    ##return e_j / (e_i + e_j)


class Sector(generic_sectors.GenericSector):
    """ Class implementing a particular sector, with attributes identifying it and methods 
    for its evaluation."""


    def __init__(self, leg_numbers, **opts):
        super(Sector, self).__init__(**opts)
        self.leg_numbers = leg_numbers


    def __call__(self, PS_point, PDGs, counterterm_index=-1, input_mapping_index=-1, sector_type=0):
        """ Given a kinematic configuration (PS_point is a *dict* here) and flavors of external states, returns the sector weight (float)
        from this sector function.
        Notice that the sectoring function will be called for both the real(s) event and *also* with the reduced
        kinematics of the counterevent, in which case it may be necessary to know which counterterm this reduced
        kinematics comes from, so that the counterterm_index is specified (not the full counterterm, because
        this function should eventually be coded at low-level and do minimal logic, so please do all the logic 
        the sector generator function. -1 for the counterterm_index implies no counterterm.
        The flavor mapping is set to None if a local counterterm is considered, and to the corresponding input
        mapping index if an integrated counterterm is considered.
        The sector_type variable determines if the standard (0) / soft (1) / collinear (2) sector function is called
        """

        q = PS_point[1] + PS_point[2]

        p_sector = [PS_point[l] for l in self.leg_numbers]

        if sector_type == 0:
            # standard sector function
            sector_weight = get_sector_wgt(q, p_sector)

            # now normalise it
            norm = 0.
            for (ii, jj) in self.all_sector_list:
                p_ii = PS_point[ii]
                p_jj = PS_point[jj]
                norm += get_sector_wgt(q, [p_ii, p_jj])

        elif sector_type == 1:
            # soft sector function
            sector_weight = get_sector_wgt_S(q, p_sector)
            # now normalise it
            norm = 0.
            for (ii, jj) in self.all_sector_list:
                # the sum runs only on the sectors with the first leg
                # equal to the one at hand
                if ii != self.leg_numbers[0]:
                    continue
                p_ii = PS_point[ii]
                p_jj = PS_point[jj]
                norm += get_sector_wgt_S(q, [p_ii, p_jj])

        elif sector_type == 2:
            # collinear sector function
            sector_weight = get_sector_wgt_C(q, p_sector)
            # it is already with the correct normalisation
            norm = 1.

        return sector_weight / norm


    def __str__(self):
        """ String representation of this sector. """

        return "(%s)"%(','.join('%d'%ln for ln in self.leg_numbers))


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

        # reeturn None if there are zero unresolved particles (virtual)

        if contrib_definition.n_unresolved_particles == 0:
            return None

        model = defining_process.get('model')
        initial_state_PDGs, final_state_PDGs = defining_process.get_cached_initial_final_pdgs()
        all_PDGs = initial_state_PDGs, final_state_PDGs

        leglist = defining_process.get('legs')

        all_sectors = []

        pert_dict = fks_common.find_pert_particles_interactions(model)
        colorlist = [model['particle_dict'][l['id']]['color'] for l in leglist]
        # the following is adapted from the very old FKS_from_real implementation
        # in the FKS formalism, i_fks is the parton associated with the soft singularities
        # however, the sector_weight function is agnostic on which parton
        # generates soft singularities
        fks_j_from_i = {}
        for i, col_i in zip(leglist, colorlist):
            if col_i == 1:
                continue
            if not i.get('state'):
                continue
            fks_j_from_i[i.get('number')] = [] # not strictly needed

            for j, col_j in zip(leglist, colorlist):
                if col_j == 1:
                    continue
                if j.get('number') == i.get('number') :
                    continue
                # if i is not a gluon, then j must not be a final state gluon
                if i['id'] != 21 and j['id'] == 21 and j['state']:
                    continue
                # if j and i are quarks and antiquark in the final state, let j be the quark
                #   this is needed in order to comply with the fct combine_ij inside fks_common
                if i['id'] == -j['id'] and j['state']:
                    if j['id'] < 0:
                        continue

                ijlist = fks_common.combine_ij(fks_common.to_fks_leg(i, model),
                                               fks_common.to_fks_leg(j, model),
                                               model, pert_dict)

                for ij in ijlist:
                    # copy the defining process, remove i and j
                    # and replace them by ij.
                    new_process = copy.copy(defining_process)
                    # this is a temporary hack waiting that squared_orders for
                    #  defining_process are correctly passed
                    ##if set(new_process['squared_orders'].values()) == set([0,]):
                    # MZ this may not be optimal, but keep for the moment
                    new_process['squared_orders'] = {}

                    new_leglist = copy.copy(leglist)
                    new_leglist[min([leglist.index(i), leglist.index(j)])] = ij
                    new_leglist.pop(max([leglist.index(i), leglist.index(j)]))
                    new_process['legs'] = new_leglist
                    if diagram_generation.Amplitude(new_process).get('diagrams'):
                        fks_j_from_i[i.get('number')].append(j.get('number'))
                        a_sector = {
                            'sector': None,
                            'counterterms': None,
                            'integrated_counterterms': None
                        }
                        a_sector['sector'] = Sector(leg_numbers=(i.get('number'), j.get('number')))
                        # keep track of the masses
                        a_sector['sector'].masses = (model.get('particle_dict')[i.get('id')]['mass'],
                                                     model.get('particle_dict')[j.get('id')]['mass'])
                        all_sectors.append(a_sector)
                        logger.info('NLO sector found, legs %d, %d' % a_sector['sector'].leg_numbers)

        if not all_sectors:
            logger.critical('WARNING, no sectors found for %s' % defining_process.nice_string())

        # up to here we have identified all the sectors.
        #  Now for each sector we need to find the corresponding counterterms
        #  We also need to add the variable all_sector_list to each of
        #  them, containing the list of all sectors (needed to normalise
        #  the weight function
        all_sector_list = [s['sector'].leg_numbers for s in all_sectors]
        all_sector_mass_list = [s['sector'].masses for s in all_sectors]

        for s in all_sectors:
            s['sector'].all_sector_list = all_sector_list
            s['sector'].all_sector_mass_list = all_sector_mass_list

            if counterterms is not None:
                s['counterterms'] = []
                for i_ct, ct in enumerate(counterterms):
                    current = ct.nodes[0].current
                    singular_structure = current.get('singular_structure').substructures[0]
                    all_legs = singular_structure.get_all_legs()
                    if singular_structure.name()=='S':
                        if all_legs[0].n == s['sector'].leg_numbers[0]: # should match to "i"
                            s['counterterms'].append(i_ct)
                    if singular_structure.name()=='C':
                        if not singular_structure.substructures:
                            # pure-collinear CT: include if the legs match those of the sector
                            if sorted([l.n for l in all_legs]) == sorted(s['sector'].leg_numbers):
                                s['counterterms'].append(i_ct)
                        else:
                            #soft-collinear CT: include only if, on top of the previous condition,
                            #  the soft leg matches the first sector leg
                            if sorted([l.n for l in all_legs]) == sorted(s['sector'].leg_numbers) and \
                               singular_structure.substructures[0].legs[0].n == s['sector'].leg_numbers[0]:
                                s['counterterms'].append(i_ct)

            # Irrelevant if this NLO example, but let me specify all of them explicitly so as to make the strucuture clear.
            if integrated_counterterms is not None:
                s['integrated_counterterms'] = {}
                for i_ct, ct in enumerate(integrated_counterterms):
                    # For now enable all integrated counterterms. Notice that the value None in this dictionary
                    # is interpreted as all input mappings contributing, but for the sake of example here
                    # we list explicitly each index.
                    s['integrated_counterterms'][i_ct] = range(len(ct['input_mappings']))

        return all_sectors
