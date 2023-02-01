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


#gl
import madgraph.interface.madgraph_interface as interface
import madgraph.iolibs.export_v4 as export
import madgraph.iolibs.file_writers as writers
import dipole_current
import torino_config
import factors_and_cuts
import recoiler_function
import colored_partons
import os
pjoin = os.path.join


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
        # print('p_sector : ' + str(self.leg_numbers))
        # print('Sector list : ' + str(self.all_sector_list))

        if sector_type == 0:
            # standard sector function
            # sigma_ij
            sector_weight = get_sector_wgt(q, p_sector)
            # add sigma_ji if j is final
            if self.leg_numbers[1] > 2:
                sector_weight += get_sector_wgt(q, [p_sector[1],p_sector[0]])
            # now normalise it
            norm = 0.
            for (ii, jj) in self.all_sector_list:
                #print('AAA - ii : ' + str(ii) + '; ' + 'jj : ' + str(jj))
                # ii runs over final state particles only
                if ii <=2:
                    raise MadEvent7Error('WARNING, sector index ii cannot be %s' % ii)
                p_ii = PS_point[ii]
                p_jj = PS_point[jj]
                # logger.info('AAA - p_ii : ' + str(ii) + '; ' + str(p_ii))
                # logger.info('AAA - p_jj : ' + str(jj) + '; ' + str(p_jj))
                if jj > 2:
                    norm += get_sector_wgt(q, [p_ii, p_jj]) + get_sector_wgt(q, [p_jj, p_ii])
                else:
                    norm += get_sector_wgt(q, [p_ii, p_jj])



        elif sector_type == 11:
            # soft sector function sigma_ij_s
            sector_weight = get_sector_wgt_S(q, p_sector)
            # now normalise it
            norm = 0.
            for (ii, jj) in self.all_sector_list:
                # print('BBB1 - ii : ' + str(ii) + '; ' + 'jj : ' + str(jj))
                # the sum runs only on the sectors with the first leg
                # equal to the one at hand
                if ii != self.leg_numbers[0] and jj != self.leg_numbers[0]:
                # if ii != self.leg_numbers[0]:
                    # print('Continuing case')
                    continue
                # print('BBB1 passed - ii : ' + str(ii) + '; ' + 'jj : ' + str(jj))
                # ii runs over final state particles only
                if ii <=2:
                    raise MadEvent7Error('WARNING, sector index ii cannot be %s' % ii)
                p_ii = PS_point[ii]
                p_jj = PS_point[jj]
                #logger.info('AAA - S p_ii : ' + str(ii) + '; ' + str(p_ii))
                # logger.info('AAA - S p_jj : ' + str(jj) + '; ' + str(p_jj))
                norm += get_sector_wgt_S(q, [p_ii, p_jj])
                # print('Norm : ' + str(norm))


        elif sector_type == 12:
            # soft sector function sigma_ji_s
            sector_weight = get_sector_wgt_S(q, [p_sector[1],p_sector[0]])
            # now normalise it
            norm = 0.
            for (ii, jj) in self.all_sector_list:
                # print('BBB2 - ii : ' + str(ii) + '; ' + 'jj : ' + str(jj))
                # the sum runs only on the sectors with the first leg
                # equal to the one at hand
                if ii != self.leg_numbers[1] and jj != self.leg_numbers[1]:
                # if jj != self.leg_numbers[1]:
                    continue
                # ii runs over final state particles only
                if ii <=2:
                    raise MadEvent7Error('WARNING, sector index ii cannot be %s' % ii)
                p_ii = PS_point[ii]
                p_jj = PS_point[jj]
                #logger.info('AAA - S p_ii : ' + str(ii) + '; ' + str(p_ii))
                #logger.info('AAA - S p_jj : ' + str(jj) + '; ' + str(p_jj))
                norm += get_sector_wgt_S(q, [p_jj, p_ii]) 


        elif sector_type == 2:
            # collinear sector function sigma_ij_c
            sector_weight = get_sector_wgt_C(q, p_sector)
            # add sigma_ji_c for j final
            if self.leg_numbers[1] > 2:
                sector_weight += get_sector_wgt_C(q, [p_sector[1],p_sector[0]])
            # now normalise it
            norm = 0.
            for (ii, jj) in self.all_sector_list:
                # the sum runs only on the sectors with the two legs
                # equal to the ones at hand
                if not (ii in self.leg_numbers and jj in self.leg_numbers):
                    continue
                # ii runs over final state particles only
                if ii <=2:
                    raise MadEvent7Error('WARNING, sector index ii cannot be %s' % ii)
                p_ii = PS_point[ii]
                p_jj = PS_point[jj]
                # logger.info('AAA - C p_ii : ' + str(ii) + '; ' + str(p_ii))
                # logger.info('AAA - C p_jj : ' + str(jj) + '; ' + str(p_jj))
                if jj <= 2:
                    norm = sector_weight
                else:
                    norm += get_sector_wgt_C(q, [p_ii, p_jj]) + get_sector_wgt_C(q, [p_jj, p_ii])


        elif sector_type == 3:
            # soft-collinear sector function
            sector_weight = get_sector_wgt_CS(q, p_sector)
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

        # return None if there are zero unresolved particles (virtual)

        if contrib_definition.n_unresolved_particles == 0:
            return None

        model = defining_process.get('model')

        initial_state_PDGs, final_state_PDGs = defining_process.get_cached_initial_final_pdgs()
        all_PDGs = initial_state_PDGs, final_state_PDGs

        leglist = defining_process.get('legs')
        #gl
        #print(defining_process.shell_string())
        #print(contrib_definition.get_shell_name())
        #print(contrib_definition.process_definition.get('id'))

        all_sectors = []
        all_sector_legs = []
        all_sector_id_legs = []

        pert_dict = fks_common.find_pert_particles_interactions(model)
        colorlist = [model['particle_dict'][l['id']]['color'] for l in leglist]

## Ez
#        logger.info("sectors.py: SectorGenerator,__call__")
#        for i,col_i in zip(leglist,colorlist):
#            logger.info("leg: "+str(i)+"\n"+str(col_i))

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
                        #print('Leg number : ' + str(a_sector['sector']))
                        #gl
                        all_sector_legs.append(i.get('number'))
                        all_sector_legs.append(j.get('number'))
                        # keep track of the masses
                        a_sector['sector'].masses = (model.get('particle_dict')[i.get('id')]['mass'],
                                                     model.get('particle_dict')[j.get('id')]['mass'])
                        #print('Masses : ' + str(a_sector['sector'].masses))
# gl
                        # keep track of the particles' identity
                        a_sector['sector'].id = (i.get('id'), j.get('id'))
                        all_sector_id_legs.append(i.get('id'))
                        all_sector_id_legs.append(j.get('id'))
                        #print('Identities : ' + str(a_sector['sector'].id))

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

        #gl
        all_local_counterterms_list = []
        necessary_ct_list = [0] * (5*len(all_sectors))
        i = 0
        for s in all_sectors:
            i = i * 5
            #print('s in sectors : ' + str(s))
            s['sector'].all_sector_list = all_sector_list
            s['sector'].all_sector_mass_list = all_sector_mass_list
            # gl
            s['sector'].all_sector_id_list = all_sector_id_list
            #print('s all_sector_list : ' + str(all_sector_list))
            #print('s all_sector_id_list : ' + str(all_sector_id_list))

            if counterterms is not None:
                s['counterterms'] = []
                for i_ct, ct in enumerate(counterterms):
                    #print('i_ct + ct : ' + str(i_ct) + ' and ' + str(ct))
                    current = ct.nodes[0].current
                    singular_structure = current.get('singular_structure').substructures[0]
                    all_legs = singular_structure.get_all_legs()
                    if singular_structure.name()=='S':
                        if all_legs[0].n == s['sector'].leg_numbers[0]: # should match to "i"
                            s['counterterms'].append(i_ct)
                            necessary_ct_list[i] = 1
# # # gl                        
                        if s['sector'].id[0] == 21 and s['sector'].id[1] == 21 and all_legs[0].n == s['sector'].leg_numbers[1]: # should match to "j"
                            s['counterterms'].append(i_ct)
                            necessary_ct_list[i+1] = 1

                    if singular_structure.name()=='C':
                        if not singular_structure.substructures:
                            # pure-collinear CT: include if the legs match those of the sector
                            if sorted([l.n for l in all_legs]) == sorted(s['sector'].leg_numbers):
                                s['counterterms'].append(i_ct)
                                necessary_ct_list[i+2] = 1
                        else:
                            #soft-collinear CT: include only if, on top of the previous condition,
                            #  the soft leg matches the first sector leg
                            if sorted([l.n for l in all_legs]) == sorted(s['sector'].leg_numbers) and \
                               singular_structure.substructures[0].legs[0].n == s['sector'].leg_numbers[0]:
                                s['counterterms'].append(i_ct)
                                necessary_ct_list[i+3] = 1
# # gl
                            if s['sector'].id[0] == 21 and s['sector'].id[1] == 21:
                                if sorted([l.n for l in all_legs]) == sorted(s['sector'].leg_numbers) and \
                                    singular_structure.substructures[0].legs[0].n == s['sector'].leg_numbers[1]:
                                    s['counterterms'].append(i_ct)
                                    necessary_ct_list[i+4] = 1

                all_local_counterterms_list.append(s['counterterms'])

            # index of necessary_ct_list    
            i = i + 1

            # Irrelevant if this NLO example, but let me specify all of them explicitly so as to make the strucuture clear.
            if integrated_counterterms is not None:
                s['integrated_counterterms'] = {}
                for i_ct, ct in enumerate(integrated_counterterms):
                    # For now enable all integrated counterterms. Notice that the value None in this dictionary
                    # is interpreted as all input mappings contributing, but for the sake of example here
                    # we list explicitly each index.
                    s['integrated_counterterms'][i_ct] = range(len(ct['input_mappings']))


######################################### Write fortran template files #############################################  

        # Set writer
        writer = writers.FortranWriter

        # TODO: point to the right process directory
        dirpath = pjoin("/Users/giovannilimatola/Desktop/Fisica/Lavori/madnklo/eejj/NLO_R_x_R_epem_guux_1")
        dirpath = pjoin(dirpath, 'SubProcesses', \
                       "P%s" % defining_process.shell_string())

        
######### Write NLO_K_isec_jsec.f

        # Set replace_dict for NLO_K_isec_jsec.f
        replace_dict_ct = {}
        for i in range(0,len(all_sector_list)):
            list_M2 = []
            isec = all_sector_list[i][0]
            jsec = all_sector_list[i][1]
            id_isec = all_sector_id_list[i][0]
            id_jsec = all_sector_id_list[i][1]
            # Check isec != jsec
            if isec == jsec:
                raise MadEvent7Error('Wrong sector indices %d,%d!' % (isec,jsec))

            replace_dict_ct['isec'] = isec
            replace_dict_ct['jsec'] = jsec
 
            if necessary_ct_list[i*5] == 1:
                if id_isec != 21:
                    raise MadEvent7Error('%d is not a gluon!' % isec)
                list_M2.append('KS=KS+M2_S(isec,xs,xp,wgt,WsumSi,xj,nitR,1d0,ierr)\n')
                list_M2.append('#\n')
            if necessary_ct_list[i*5+1] == 1:
                if id_jsec != 21:
                    raise MadEvent7Error('%d is not a gluon!' % jsec)
                list_M2.append('KS=KS+M2_S(jsec,xs,xp,wgt,WsumSj,xj,nitR,1d0,ierr)\n')
                list_M2.append('#\n')
            if necessary_ct_list[i*5+2] == 1:
                # Loop over sectors with final state particles only
                if isec > 2 and jsec > 2:
                    # Extract the reference particle leg from recoiler_function.py
                    iref_leg = recoiler_function.get_recoiler(defining_process,(isec,jsec))
                    iref = iref_leg.get('number')
                    replace_dict_ct['iref'] = iref
                    # Check irec validity
                    if (isec == iref) or (jsec == iref):
                        raise MadEvent7Error('Wrong recoiler %d,%d,%d!' % (isec,jsec,iref))
                    # Write an identified M2_H_C_F*F* for each (**) flavour couple 
                    if id_isec == 21 and id_jsec == 21:
                        list_M2.append('KHC=KHC+M2_H_C_FgFg(isec,jsec,%d,xs,xp,xsb,xpb,wgt,xj,nitR,1d0,ierr)' % iref)
                    elif id_isec == 21 and id_jsec != 21: # if there is a gluon in sector, it is always in the first position
                        list_M2.append('KHC=KHC+M2_H_C_FgFq(isec,jsec,%d,xs,xp,xsb,xpb,wgt,xj,nitR,1d0,ierr)' % iref)
                    else:
                        list_M2.append('KHC=KHC+M2_H_C_FqFqx(isec,jsec,%d,xs,xp,xsb,xpb,wgt,xj,nitR,1d0,ierr)' % iref)
                # Loop over sectors with at least one initial state particle
                if isec <= 2 or jsec <= 2:
                    continue

            # soft-collinear kernels
            #if str_cts[i*10+6] == 1:
            #    str_M2.append('')
            #if str_cts[i*10+8] == 1:
            #    str_M2.append('')

                str_M2 = " ".join(list_M2)
                replace_dict_ct['str_M2'] = str_M2

            filename = pjoin(dirpath, 'NLO_K_%d_%d.f' % (isec, jsec))
            file = open("/Users/giovannilimatola/Desktop/Fisica/Lavori/madnklo/tmp_fortran/tmp_files/NLO_K_template.f").read()
            file = file % replace_dict_ct
            writer(filename).writelines(file)

######## Write masses.inc ########
        replace_dict={}
        scale = model.get('parameter_dict')['MU_R']
        mz = model.get('parameter_dict')['mdl_MZ']
        mt = model.get('parameter_dict')['mdl_MT']
        mc = model.get('parameter_dict')['mdl_MC']
        mb = model.get('parameter_dict')['mdl_MB']
        me = model.get('parameter_dict')['mdl_Me']
        mmu = model.get('parameter_dict')['mdl_MM']
        mta = model.get('parameter_dict')['mdl_MTA']
        replace_dict['MU_R'] = scale
        replace_dict['mdl_MZ'] = mz
        replace_dict['mdl_MT'] = mt
        replace_dict['mdl_MC'] = mc
        replace_dict['mdl_MB'] = mb
        replace_dict['mdl_Me'] = me
        replace_dict['mdl_MM'] = mmu
        replace_dict['mdl_MTA'] = mta
        
        file = """ \
          double precision mur
          double precision me,mmu,mta
          double precision mc,mt,mb
          double precision mz,mw
          parameter (mz = %(mdl_MZ)lfd0)
          parameter (me = %(mdl_Me)lfd0, mmu = %(mdl_MM)lfd0, mta = %(mdl_MTA)lfd0)
          parameter (mc = %(mdl_MC)lfd0, mt = %(mdl_MT)lfd0, mb = %(mdl_MB)lfd0)
          parameter (mur = %(MU_R)lfd0)
          common/model/me,mmu,mta,mc,mt,mb,mz,mw""" % replace_dict

        filename = pjoin(dirpath, 'masses_factors.inc')
        writer(filename).writelines(file)

        
######### Write damping_factors.inc

        # Set replace_dict 
        replace_dict = {}
        dfactors = factors_and_cuts.damping_factors
        replace_dict['alpha'] = dfactors[0]
        replace_dict['beta_FF'] = dfactors[1]
        replace_dict['beta_FI'] = dfactors[2]
        replace_dict['beta_IF'] = dfactors[3]
        replace_dict['beta_II'] = dfactors[4]

        file = """ \
          double precision alpha, beta_FF, beta_FI, beta_IF, beta_II
          parameter (alpha = %(alpha)dd0)
          parameter (beta_FF = %(beta_FF)dd0)
          parameter (beta_FI = %(beta_FI)dd0)
          parameter (beta_IF = %(beta_IF)dd0)
          parameter (beta_II = %(beta_II)dd0)""" % replace_dict

        filename = pjoin(dirpath, 'damping_factors.inc')
        writer(filename).writelines(file)


######### Write colored_partons.inc

        # Set replace_dict 
        replace_dict = {}
        colored_legs = list(colored_partons.get_colored_legs(defining_process).values())
        list_colored_legs = []
        isNLOQCDparton = ['.false.'] * (len(leglist))
        for i in range(0,len(colored_legs)):
            list_colored_legs.append(colored_legs[i][0])
            isNLOQCDparton[colored_legs[i][0]-1] = '.true.'
        replace_dict['len_leglist'] = len(leglist)
        replace_dict['isNLOQCDparton'] = str(isNLOQCDparton).replace('[','').replace(']','').replace(' ','').replace("'","")

        file = """ \
          logical isNLOQCDparton(%(len_leglist)d)
          data isNLOQCDparton/%(isNLOQCDparton)s/""" % replace_dict

        filename = pjoin(dirpath, 'colored_partons.inc')
        writer(filename).writelines(file)






#        # Set replace_dict for sector.inc
#        replace_dict = {}
#        replace_dict['singular_sec'] = len(all_sector_list)
#        replace_dict['all_sector_legs'] = str(all_sector_legs).replace('[','').replace(']','').replace(' ','')
#        replace_dict['all_sector_id_legs'] = str(all_sector_id_legs).replace('[','').replace(']','').replace(' ','')
#        replace_dict['necessary_ct_list'] = str(necessary_ct_list).replace('[','').replace(']','').replace(' ','')

#        filename = pjoin(dirpath, 'sector_test.inc')
#        file = open("/home/gloria/Desktop/madnklo/tmp_fortran/template_files/sector.inc").read()
#        file = file % replace_dict
#        writer(filename).writelines(file)

        #gl
        """# Wrote sector_template file
        writer = writers.FortranWriter
        dirpath = "/Users/giovannilimatola/Desktop/Fisica/Lavori/madnklo/tmp_fortran/template_files"
        filename = pjoin(dirpath, 'sector.inc')

        replace_dict = {}
        replace_dict['all_sector_list'] = all_sector_list
        replace_dict['all_sector_id_list'] = all_sector_id_list
        replace_dict['all_local_counterterms_list'] = all_local_counterterms_list

        file = open("/Users/giovannilimatola/Desktop/Fisica/Lavori/madnklo/tmp_fortran/template_files/sector_template.inc").read()
        file = file % replace_dict
        writer(filename).writelines(file)"""


        return all_sectors
