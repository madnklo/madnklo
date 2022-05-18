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



class MadEvent7Error(Exception):
    pass


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
    e_j = s_qj / s
    w_ij = s * s_ij / s_qi / s_qj
    # logger.info('s : ' + str(s) + '; ' + 'q : ' + str(q))
    # logger.info('BBB - s_ij : ' + str(s_ij))
    # logger.info('BBB - s_qi : ' + str(s_qi))
    # logger.info('BBB - s_qj : ' + str(s_qj))
    # logger.info('BBB - e_i : ' + str(e_i))
    # logger.info('BBB - w_ij : ' + str(w_ij))

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
    # logger.info('BBB - S s_ij : ' + str(s_ij))
    # logger.info('BBB - S s_qi : ' + str(s_qi))
    # logger.info('BBB - S s_qj : ' + str(s_qj))
    # logger.info('BBB - S w_ij : ' + str(w_ij))
    return 1 / w_ij


def get_sector_wgt_C(q, p_sector):
    """ the collinear sector function sigma_ij of the Torino paper (eq. 2.15right)
     without the normalisation.
     - q is the total momentum of the incoming particles
     - p_sector is a list of the momenta of the particles defining the sector
    """
    s = q.square()
    # s_ij = 2 * p_sector[0].dot(p_sector[1])
    # s_qj = 2 * p_sector[1].dot(q)
    s_qi = 2 * p_sector[0].dot(q)
    e_i = s_qi / s
    # e_j = s_qj / s
    # logger.info('BBB - C s_ij : ' + str(s_ij))
    # logger.info('BBB - C s_qj : ' + str(s_qj))
    # logger.info('BBB - C e_j : ' + str(e_j))
    # logger.info('BBB - C s_qi : ' + str(s_qi))
    # logger.info('BBB - C e_i : ' + str(e_i))
    return 1. / e_i


def get_sector_wgt_CS(q, p_sector):
    """ the soft-collinear sector function sigma_ij of the Torino paper (eq. 2.15right)
     without the normalisation.
     - q is the total momentum of the incoming particles
     - p_sector is a list of the momenta of the particles defining the sector
     This function just returns 1
    """
    return 1


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
        The sector_type variable determines if the standard (0) / soft (1) / collinear (2) / soft-collinear(3)
        sector function is called
        """

        q = PS_point[1] + PS_point[2]

        p_sector = [PS_point[l] for l in self.leg_numbers]
        p_sector_ids = self.all_sector_id_list
        # print('p_sector identities : ' + str(self.all_sector_id_list))


        if sector_type == 0:
            # standard sector function
            # sigma_ij
            sector_weight = get_sector_wgt(q, p_sector)
            # add sigma_ji
            sector_weight += get_sector_wgt(q, [p_sector[1],p_sector[0]])
            # now normalise it
            norm = 0.
            for (ii, jj) in self.all_sector_list:
                # print('AAA - ii : ' + str(ii) + '; ' + 'jj : ' + str(jj))
                # ii runs over final state particles only
                if ii <=2:
                    raise MadEvent7Error('WARNING, sector index ii cannot be %s' % ii)
                p_ii = PS_point[ii]
                p_jj = PS_point[jj]
                # logger.info('AAA - p_ii : ' + str(ii) + '; ' + str(p_ii))
                # logger.info('AAA - p_jj : ' + str(jj) + '; ' + str(p_jj))
                # logger.info('AAA - get_sector_weight : ' + str(get_sector_wgt(q, [p_ii, p_jj])))
                #print('Sector functions : ' + str(get_sector_wgt(q, [p_ii, p_jj]) + get_sector_wgt(q, [p_jj, p_ii])))
                norm += get_sector_wgt(q, [p_ii, p_jj]) + get_sector_wgt(q, [p_jj, p_ii])
                #print('Norm : ' + str(norm))
            # print('AAA - sector_weight : ' + str(sector_weight))
            # print('AAA - norm : ' + str(norm))
            # logger.info('AAA - norm : ' + str(norm))
            # logger.info('AAA - sector_weight : ' + str(sector_weight))
            # logger.info('AAA - norm : ' + str(norm))

        elif sector_type == 11:
            # soft sector function sigma_ij_s
            #sector_weight = get_sector_wgt(q, p_sector)
            sector_weight = get_sector_wgt_S(q, p_sector)
            # now normalise it
            norm = 0.
            # norm_2 = 0.
            for (ii, jj) in self.all_sector_list:
                # print('BBB1 - ii : ' + str(ii) + '; ' + 'jj : ' + str(jj))
                # the sum runs only on the sectors with the first leg
                # equal to the one at hand
                #print('leg numbers : ' + str(self.leg_numbers[1].get('id')))
                if ii != self.leg_numbers[0] and jj != self.leg_numbers[0]:
                    continue
                # ii runs over final state particles only
                if ii <=2:
                    raise MadEvent7Error('WARNING, sector index ii cannot be %s' % ii)
                p_ii = PS_point[ii]
                p_jj = PS_point[jj]
                #logger.info('AAA - S p_ii : ' + str(ii) + '; ' + str(p_ii))
                # logger.info('AAA - S p_jj : ' + str(jj) + '; ' + str(p_jj))
                # logger.info('BBB1 - S get_sector_weight_S : ' + str(get_sector_wgt_S(q, [p_ii, p_jj])))
                # norm += get_sector_wgt(q, [p_ii, p_jj]) + get_sector_wgt(q, [p_jj, p_ii])
                norm += get_sector_wgt_S(q, [p_ii, p_jj])


        elif sector_type == 12:
            # soft sector function sigma_ji_s
            # sector_weight = get_sector_wgt(q, [p_sector[1],p_sector[0]])
            sector_weight = get_sector_wgt_S(q, [p_sector[1],p_sector[0]])
            # now normalise it
            norm = 0.
            for (ii, jj) in self.all_sector_list:
                # print('BBB2 - ii : ' + str(ii) + '; ' + 'jj : ' + str(jj))
                # the sum runs only on the sectors with the first leg
                # equal to the one at hand
                #print('leg numbers : ' + str(self.leg_numbers[1].get('id')))
                if ii != self.leg_numbers[1] and jj != self.leg_numbers[1]:
                    continue
                # ii runs over final state particles only
                if ii <=2:
                    raise MadEvent7Error('WARNING, sector index ii cannot be %s' % ii)
                p_ii = PS_point[ii]
                p_jj = PS_point[jj]
                #logger.info('AAA - S p_ii : ' + str(ii) + '; ' + str(p_ii))
                # logger.info('AAA - S p_jj : ' + str(jj) + '; ' + str(p_jj))
                #logger.info('BBB2 - S get_sector_weight_S : ' + str(get_sector_wgt_S(q, [p_ii, p_jj])))
                # norm += get_sector_wgt(q, [p_jj, p_ii]) + get_sector_wgt(q, [p_ii, p_jj])
                norm += get_sector_wgt_S(q, [p_jj, p_ii]) 


        elif sector_type == 2:
            # collinear sector function sigma_ij_c
            sector_weight = get_sector_wgt_C(q, p_sector)
            # add sigma_ji_c
            sector_weight += get_sector_wgt_C(q, [p_sector[1],p_sector[0]])

            # now normalise it
            norm = 0.
            for (ii, jj) in self.all_sector_list:
                # the sum runs only on the sectors with the two legs
                # equal to the ones at hand
                # logger.info('AAA - C self.all_sector_list : ' + str(self.all_sector_list))
                # logger.info('AAA - C self.leg_numbers : ' + str(self.leg_numbers))
                # logger.info('AAA - jj : ' + str(jj))
                if not (ii in self.leg_numbers and jj in self.leg_numbers):
                    continue
                # ii runs over final state particles only
                if ii <=2:
                    raise MadEvent7Error('WARNING, sector index ii cannot be %s' % ii)
                p_ii = PS_point[ii]
                p_jj = PS_point[jj]
                # logger.info('AAA - C p_ii : ' + str(ii) + '; ' + str(p_ii))
                # logger.info('AAA - C p_jj : ' + str(jj) + '; ' + str(p_jj))
                # norm += get_sector_wgt_C(q, [p_ii, p_jj])
                # logger.info('AAA - C norm dentro circuito : ' + str(norm) )
                if jj <= 2:
                    #logger.info('AAA - C jj : ' + str(jj))
                    norm = sector_weight
                else:
                    norm += get_sector_wgt_C(q, [p_ii, p_jj]) + get_sector_wgt_C(q, [p_jj, p_ii])
            # logger.info('AAA - C sector_weight : ' + str(sector_weight))
            # logger.info('AAA - C norm : ' + str(norm))
            # logger.info('AAA - C result : ' + str(sector_weight/norm))

        elif sector_type == 3:
            # soft-collinear sector function
            sector_weight = get_sector_wgt_CS(q, p_sector)
            # it is already with the correct normalisation
            norm = 1.

            # logger.info('AAA - SC sector_weight : ' + str(sector_weight))
            # logger.info('AAA - SC norm : ' + str(norm))
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

        # return None if there are zero unresolved particles (virtual)

        if contrib_definition.n_unresolved_particles == 0:
            return None

        model = defining_process.get('model')
        initial_state_PDGs, final_state_PDGs = defining_process.get_cached_initial_final_pdgs()
        all_PDGs = initial_state_PDGs, final_state_PDGs

        leglist = defining_process.get('legs')
        # print('leglist : ' + str(leglist))

        all_sectors = []

        pert_dict = fks_common.find_pert_particles_interactions(model)
        colorlist = [model['particle_dict'][l['id']]['color'] for l in leglist]

# Ez
        logger.info("sectors.py: SectorGenerator,__call__")
        for i,col_i in zip(leglist,colorlist):
            logger.info("leg: "+str(i)+"\n"+str(col_i))

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
# gl
                # if both i and j are gluons, then keep just the case in which i (number) < j (number)
                if i['id'] == 21 and j['id'] == 21 and j['state']:
                    if j.get('number') < i.get('number') :
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
                        print('Leg number : ' + str(a_sector['sector']))
                        # keep track of the masses
                        a_sector['sector'].masses = (model.get('particle_dict')[i.get('id')]['mass'],
                                                     model.get('particle_dict')[j.get('id')]['mass'])
                        print('Masses : ' + str(a_sector['sector'].masses))
# gl
                        # keep track of the particles' identity
                        a_sector['sector'].id = (i.get('id'), j.get('id'))
                        print('Identities : ' + str(a_sector['sector'].id))

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
        # gl
        all_sector_id_list = [s['sector'].id for s in all_sectors]

        for s in all_sectors:
            s['sector'].all_sector_list = all_sector_list
            s['sector'].all_sector_mass_list = all_sector_mass_list
            # gl
            s['sector'].all_sector_id_list = all_sector_id_list

            if counterterms is not None:
                s['counterterms'] = []
                for i_ct, ct in enumerate(counterterms):
                    current = ct.nodes[0].current
                    singular_structure = current.get('singular_structure').substructures[0]
                    all_legs = singular_structure.get_all_legs()
                    if singular_structure.name()=='S':
                        if all_legs[0].n == s['sector'].leg_numbers[0]: # should match to "i"
                            s['counterterms'].append(i_ct)
# # gl
                        if all_legs[0].n == s['sector'].leg_numbers[1]: # should match to "i"
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
# gl
                            if sorted([l.n for l in all_legs]) == sorted(s['sector'].leg_numbers) and \
                               singular_structure.substructures[0].legs[0].n == s['sector'].leg_numbers[1]:
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
