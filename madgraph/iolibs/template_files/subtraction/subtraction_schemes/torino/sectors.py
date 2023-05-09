#
# Generic classes for the sector implementation in subtraction schemes
#
import copy
import sys


import commons.generic_sectors as generic_sectors
import madgraph.various.misc as misc
import madgraph.core.diagram_generation as diagram_generation
import madgraph.fks.fks_common as fks_common
import madgraph.integrator.vectors as vectors
import logging


#gl
import madgraph.interface.madgraph_interface as interface
import madgraph.iolibs.file_writers as writers
import madgraph.core.contributions as contributions
import dipole_current
import factors_and_cuts
import recoiler_function
import colored_partons
import os
import glob
pjoin = os.path.join
from madgraph.iolibs.files import cp, ln, mv

import madgraph.iolibs.export_v4 as export_v4
import madgraph.iolibs.export_ME7 as export_ME7
import madgraph.interface.madgraph_interface as interface


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
                print('BBB1 - ii : ' + str(ii) + '; ' + 'jj : ' + str(jj))
                # the sum runs only on the sectors with the first leg
                # equal to the one at hand
                if ii != self.leg_numbers[0] and jj != self.leg_numbers[0]:
                # if ii != self.leg_numbers[0]:
                    # print('Continuing case')
                    continue
                print('BBB1 passed - ii : ' + str(ii) + '; ' + 'jj : ' + str(jj))
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
                'integrated_counterterms' : [list of tuples (counterterm indices, input_mapping_index) for integrated CT to be considered],
                'recoiler' : [recoiler leg]
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
        all_sector_recoilers = []

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
                            'integrated_counterterms': None,
                            'recoiler' : None
                        }
                        a_sector['sector'] = Sector(leg_numbers=(i.get('number'), j.get('number')))
                        a_sector['recoiler'] = recoiler_function.get_recoiler(defining_process,(i.get('number'),j.get('number')))
                        all_sector_recoilers.append(a_sector['recoiler'].get('number'))
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
        necessary_ct = [0] * (5*len(all_sectors))
        i = 0
        for s in all_sectors:
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
                            necessary_ct[i] = ct
# # # gl                        
                        if s['sector'].id[0] == 21 and s['sector'].id[1] == 21 and all_legs[0].n == s['sector'].leg_numbers[1]: # should match to "j"
                            s['counterterms'].append(i_ct)
                            necessary_ct_list[i+1] = 1
                            necessary_ct[i+1] = ct

                    if singular_structure.name()=='C':
                        if not singular_structure.substructures:
                            # pure-collinear CT: include if the legs match those of the sector
                            if sorted([l.n for l in all_legs]) == sorted(s['sector'].leg_numbers):
                                s['counterterms'].append(i_ct)
                                necessary_ct_list[i+2] = 1
                                necessary_ct[i+2] = ct
                        else:
                            #soft-collinear CT: include only if, on top of the previous condition,
                            #  the soft leg matches the first sector leg
                            if sorted([l.n for l in all_legs]) == sorted(s['sector'].leg_numbers) and \
                               singular_structure.substructures[0].legs[0].n == s['sector'].leg_numbers[0]:
                                s['counterterms'].append(i_ct)
                                necessary_ct_list[i+3] = 1
                                necessary_ct[i+3] = ct
# # gl
                            if s['sector'].id[0] == 21 and s['sector'].id[1] == 21:
                                if sorted([l.n for l in all_legs]) == sorted(s['sector'].leg_numbers) and \
                                    singular_structure.substructures[0].legs[0].n == s['sector'].leg_numbers[1]:
                                    s['counterterms'].append(i_ct)
                                    necessary_ct_list[i+4] = 1
                                    necessary_ct[i+4] = ct

                all_local_counterterms_list.append(s['counterterms'])

            # index of necessary_ct_list    
            i += 5 

            # Irrelevant if this NLO example, but let me specify all of them explicitly so as to make the strucuture clear.
            if integrated_counterterms is not None:
                s['integrated_counterterms'] = {}
                for i_ct, ct in enumerate(integrated_counterterms):
                    # For now enable all integrated counterterms. Notice that the value None in this dictionary
                    # is interpreted as all input mappings contributing, but for the sake of example here
                    # we list explicitly each index.
                    s['integrated_counterterms'][i_ct] = range(len(ct['input_mappings']))

######################################### Write fortran template files for n+1 body #############################################  

        # Set writer
        writer = writers.FortranWriter

        # TODO: point to the right process directory
        dirmadnklo=os.getcwd()
        #file = glob.glob("eejj/NLO_R_x_R_*")
        #file = glob.glob("%s/NLO_R_x_R_*" % interface.user_dir_name[0])
        #print(file)
        dirpath = pjoin(dirmadnklo,glob.glob("%s/NLO_R_x_R_*" % interface.user_dir_name[0])[0])
        #dirpath = pjoin(dirmadnklo,"eejj/NLO_R_x_R_epem_guux_1")
        dirpath_subprocesses = pjoin(dirpath, 'SubProcesses')
        dirpath = pjoin(dirpath, 'SubProcesses', \
                       "P%s" % defining_process.shell_string())


######### Import Born-level PDGs

        tmp_dirpath = pjoin(dirmadnklo,"%s/SubProcesses" % interface.user_dir_name[0])
        sys.path.append(tmp_dirpath)
        import Born_PDGs as PDGs_from_Born

        v_rec = recoiler_function.get_virtual_recoiler(PDGs_from_Born.leg_PDGs_epem_ddx)
        #print(v_rec)


######### Write all_sector_list.inc

        replace_dict = {}
        replace_dict['len_sec_list'] = len(all_sector_list)
        replace_dict['all_sector_list'] = str(all_sector_list).replace('[','').replace(']','').replace(' ','').replace('(','').replace(')','')

        file = """ \
          integer, parameter :: lensectors = %(len_sec_list)d
          integer all_sector_list(2,lensectors)
          data all_sector_list/%(all_sector_list)s/""" % replace_dict

        filename = pjoin(dirpath, 'all_sector_list.inc')
        writer(filename).writelines(file)
        
######### Write NLO_K_isec_jsec.f, NLO_Rsub_isec_jsec.f

        # Set replace_dict for NLO_K_isec_jsec.f
        replace_dict_ct = {}
        replace_dict_int_real = {}
        list_virtual = []
        for i in range(0,len(all_sector_list)):
            list_M2 = []
            list_str_defHC = []
            list_int_real = []
            mapping = []
            mapping_str = ''
            isec = all_sector_list[i][0]
            jsec = all_sector_list[i][1]
            id_isec = all_sector_id_list[i][0]
            id_jsec = all_sector_id_list[i][1]
            # Check isec != jsec
            if isec == jsec:
                raise MadEvent7Error('Wrong sector indices %d,%d!' % (isec,jsec))
            replace_dict_ct['isec'] = isec
            replace_dict_ct['jsec'] = jsec
            replace_dict_int_real['isec'] = isec
            replace_dict_int_real['jsec'] = jsec
            # iref for phase_space_npo in NLO_Rsub
            replace_dict_int_real['iref'] = all_sector_recoilers[i]
            list_virtual.append(isec)
            list_virtual.append(jsec)
            list_virtual.append(all_sector_recoilers[i])
            
            if necessary_ct_list[i*5] == 1:
                if id_isec != 21:
                    raise MadEvent7Error('%d is not a gluon!' % isec)
                list_M2.append('KS=KS+M2_S(isec,xs,xp,wgt,WsumSi,xj,nitR,1d0,ierr)\n')
                list_M2.append('if(ierr.eq.1)goto 999\n')
                list_int_real.append('# call sector function ZsumSi\n')
                list_int_real.append('call get_Z_NLO(sNLO,sCM,alpha,isec,jsec,ZsumSi,"S",ierr)\n')
                list_int_real.append('if(ierr.eq.1)goto 999\n')
            if necessary_ct_list[i*5+1] == 1:
                if id_jsec != 21:
                    raise MadEvent7Error('%d is not a gluon!' % jsec)
                list_M2.append('KS=KS+M2_S(jsec,xs,xp,wgt,WsumSj,xj,nitR,1d0,ierr)\n')
                list_M2.append('if(ierr.eq.1)goto 999\n')
                list_int_real.append('# call sector function ZsumSj\n')
                list_int_real.append('call get_Z_NLO(sNLO,sCM,alpha,jsec,isec,ZsumSj,"S",ierr)\n')
                list_int_real.append('if(ierr.eq.1)goto 999\n')
            if necessary_ct_list[i*5+2] == 1:
                # Loop over sectors with final state particles only
                if isec > 2 and jsec > 2:
                    # Extract the reference particle leg from recoiler_function.py
                    #iref_leg = recoiler_function.get_recoiler(defining_process,(isec,jsec))
                    #iref = iref_leg.get('number')
                    iref = all_sector_recoilers[i]
                    replace_dict_ct['iref'] = iref
                    # Check irec validity
                    if (isec == iref) or (jsec == iref):
                        raise MadEvent7Error('Wrong recoiler %d,%d,%d!' % (isec,jsec,iref))
                    # Write an identified M2_H_C_F*F* for each (**) flavour couple 
                    if id_isec == 21 and id_jsec == 21:
                        list_M2.append('KHC=KHC+M2_H_C_FgFg(isec,jsec,iref,xs,xp,xsb,xpb,wgt,xj,nitR,1d0,ierr)\n')
                        list_str_defHC.append('DOUBLE PRECISION M2_H_C_FgFg')
                    elif id_isec == 21 and id_jsec != 21: # if there is a gluon in sector, it is always in the first position
                        list_M2.append('KHC=KHC+M2_H_C_FgFq(isec,jsec,iref,xs,xp,xsb,xpb,wgt,xj,nitR,1d0,ierr)\n')
                        list_str_defHC.append('DOUBLE PRECISION M2_H_C_FgFq')
                    else:
                        list_M2.append('KHC=KHC+M2_H_C_FqFqx(isec,jsec,iref,xs,xp,xsb,xpb,wgt,xj,nitR,1d0,ierr)\n')
                        list_str_defHC.append('DOUBLE PRECISION M2_H_C_FqFqx')
                    list_M2.append('if(ierr.eq.1)goto 999\n')
                    #list_int_real.append('# Symmetrised sector function for collinear kernel is equal to 1\n')

                    # default mapping for final-state collinear kernels
                    # (abc) == (ijr)
                    mapping = ['ISEC', 'JSEC', 'IREF'] 
                # Loop over sectors with at least one initial state particle
                if isec <= 2 or jsec <= 2:
                    continue

            # soft-collinear kernels
            #if str_cts[i*10+6] == 1:
            #    str_M2.append('')
            #if str_cts[i*10+8] == 1:
            #    str_M2.append('')

                #specify (abc) mapping choice
                mapping_str = """ \
                    iU = %s
                    iS = %s
                    iB = %s
                    iA = 1 ! default azimuth for NLO""" % (mapping[0],mapping[1],mapping[2])

                str_defHC = " ".join(list_str_defHC)
                str_M2 = " ".join(list_M2)
                str_int_real = " ".join(list_int_real)
                replace_dict_int_real['mapping_str'] = mapping_str
                replace_dict_ct['str_defHC'] = str_defHC
                replace_dict_ct['str_M2'] = str_M2
                replace_dict_int_real['str_int_real'] = str_int_real        
                replace_dict_int_real['mapping_str'] = mapping_str       
                replace_dict_int_real['NLO_proc_str'] = str(defining_process.shell_string(schannel=True, 
                                        forbid=True, main=False, pdg_order=False, print_id = False) + '_')

            # write NLO_K
            filename = pjoin(dirpath, 'NLO_K_%d_%d.f' % (isec, jsec))
            file = open(pjoin(dirmadnklo,"tmp_fortran/tmp_files/NLO_K_template.f")).read()
            file = file % replace_dict_ct
            writer(filename).writelines(file)
            # write NLO_Rsub
            filename_int_real = pjoin(dirpath, 'NLO_Rsub_%d_%d.f' % (isec, jsec))
            file_int_real = open(pjoin(dirmadnklo,"tmp_fortran/tmp_files/NLO_Rsub_template.f")).read()
            file_int_real = file_int_real % replace_dict_int_real
            writer(filename_int_real).writelines(file_int_real)
        # # write virtual_recoilers.inc
        # str_virtual = " ".join(str(list_virtual))
        # replace_dict['len_sec_list'] = len(all_sector_list)
        # replace_dict['ijr_set'] = str_virtual.replace('[','').replace(']','').replace(' ','').replace('(','').replace(')','')
        # file = """ \
        #   integer, parameter :: lensectors = %(len_sec_list)d
        #   integer sector_particles_ijr(3,lensectors)
        #   data sector_particles_ijr/%(ijr_set)s/""" % replace_dict
        # filename = pjoin(dirpath, 'virtual_recoilers.inc')
        # writer(filename).writelines(file)


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
          double precision alpha
          double precision beta_FF,beta_FI
          double precision beta_IF,beta_II
          parameter (alpha = %(alpha)fd0)
          parameter (beta_FF = %(beta_FF)dd0)
          parameter (beta_FI = %(beta_FI)dd0)
          parameter (beta_IF = %(beta_IF)dd0)
          parameter (beta_II = %(beta_II)dd0)""" % replace_dict

        filename = pjoin(dirpath, 'damping_factors.inc')
        writer(filename).writelines(file)


######### Write colored_partons.inc

        # # Set replace_dict 
        # replace_dict = {}
        # colored_legs = list(colored_partons.get_colored_legs(defining_process).values())
        # list_colored_legs = []
        # isNLOQCDparton = ['.false.'] * (len(leglist))
        # for i in range(0,len(colored_legs)):
        #     list_colored_legs.append(colored_legs[i][0])
        #     isNLOQCDparton[colored_legs[i][0]-1] = '.true.'
        # replace_dict['len_leglist'] = len(leglist)
        # replace_dict['isNLOQCDparton'] = str(isNLOQCDparton).replace('[','').replace(']','').replace(' ','').replace("'","")

        # file = """ \
        #   logical isNLOQCDparton(%(len_leglist)d)
        #   data isNLOQCDparton/%(isNLOQCDparton)s/""" % replace_dict

        # filename = pjoin(dirpath, 'colored_partons.inc')
        # writer(filename).writelines(file)


# ######### Write leg_PDGs.inc

        # # Set replace_dict 
        # replace_dict = {}
        # leg_PDGs = []
        # leg_PDGs.append(all_PDGs[0][0])
        # leg_PDGs.append(all_PDGs[0][1])
        # for i in range(0,len(final_state_PDGs)):
        #     leg_PDGs.append(all_PDGs[1][i])

        # replace_dict['len_legPDGs'] = len(leg_PDGs)
        # replace_dict['leg_PDGs'] = str(leg_PDGs).replace('[','').replace(']','').replace(' ','').replace("'","")

        # file = """ \
        #   integer leg_PDGs(%(len_legPDGs)d)
        #   data leg_PDGs/%(leg_PDGs)s/""" % replace_dict

        # filename = pjoin(dirpath, 'leg_PDGs.inc')
        # writer(filename).writelines(file)

########################################################################

        #file_LO = glob.glob("%s/LO_*" % interface.user_dir_name[0])
        #dirpathLO = pjoin(dirmadnklo,glob.glob("%s/LO_*" % interface.user_dir_name[0])[0])
        # dirpathLO = pjoin(dirpathLO, 'SubProcesses', \
        #                "P%s" % counterterms[0].current.shell_string())

        # # Extract for each counterterm the respecitve underlying Born files (matrix_*.f, spin_correlations.inc) 
        # Born_processes = []
        # for i, ct in enumerate(counterterms):
        #     new_proc = ct.current.shell_string(
        #         schannel=False, forbid=False, main=False, pdg_order=False, print_id = False)
        #     if new_proc not in Born_processes:
        #         Born_processes.append(new_proc)
        
        # print(Born_processes)

 
        # os.symlink( "%s/matrix_%s.f" % (dirpathLO, counterterms[0].current.shell_string(
        #     schannel=False, forbid=False, main=False, pdg_order=False, print_id = False)), 
        #     "%s/matrix_%s.f" % (dirpath, counterterms[0].current.shell_string(
        #     schannel=False, forbid=False, main=False, pdg_order=False, print_id = False)) )
        # os.symlink(dirpathLO + '/spin_correlations.inc' , dirpath + '/spin_correlations.inc' )

        # # #gl
        # for i, ct in enumerate(counterterms):
        #     # return ((S(3),),)
        #     print(ct)
        #     # return "1_epem_uux"
        #     print(ct.current.shell_string())
        #     # return "epem_uux"
        #     print(ct.current.shell_string(
        #     schannel=False, forbid=False, main=False, pdg_order=False, print_id = False))
        #     # return ((-11, 11), (2, -2))
        #     print(ct.get_reduced_flavors())

######### Write NLO_IR_limits_isec_jsec.f and import underlying Born MEs and spin_correlations.inc

        dirpathLO = pjoin(dirmadnklo,glob.glob("%s/LO_*" % interface.user_dir_name[0])[0])

        replace_dict_limits = {}
        # List of necessary underlying Born strings and particle PDGs
        Born_processes = []
        Born_PDGs = []
        # List of dirpathLO of the necessary underlying Born
        path_Born_processes = []
        # Link LO files to each real process directory
        dirpathLO_head = pjoin(dirmadnklo,glob.glob("%s/LO_*" % interface.user_dir_name[0])[0])
        for i in range(0,len(all_sector_list)):
            tmp_Born_PDGs = []
            isec = all_sector_list[i][0]
            jsec = all_sector_list[i][1]

            replace_dict_limits['isec'] = isec
            replace_dict_limits['jsec'] = jsec

            replace_dict_limits['proc_prefix_S'] = 'dummy'
            replace_dict_limits['proc_prefix_H_C_FgFg'] = 'dummy'
            replace_dict_limits['proc_prefix_H_C_FgFq'] = 'dummy'
            replace_dict_limits['proc_prefix_H_C_FqFqx'] = 'dummy'

            if necessary_ct_list[i*5] == 1 or necessary_ct_list[i*5+1] == 1:
                # list of proc str permutations 'epem_ddx' for template
                uB_proc = necessary_ct[i*5].current.shell_string_user(
                            schannel=True, forbid=True, main=False, pdg_order=False, print_id = False)
                # list of proc str permutations '1_epem_ddx' for directory
                uB_proc_str_1 = necessary_ct[i*5].current.shell_string_user()
                for j in range(0,len(uB_proc)):
                    dirpathLO = pjoin(dirpathLO_head, 'SubProcesses', "P%s" % uB_proc_str_1[j])
                    if os.path.exists(dirpathLO):
                        replace_dict_limits['proc_prefix_S'] = uB_proc[j]
                        tmp_Born_PDGs.append((uB_proc[j],isec,jsec))
                        if uB_proc[j] not in Born_processes:
                            Born_processes.append(uB_proc[j])
                            path_Born_processes.append(dirpathLO)
                        break
                    # grouped subprocesses have no a specific LO directory
                    if j == len(uB_proc) - 1:
                        extra_uB_proc = uB_proc[0]
                        replace_dict_limits['proc_prefix_S'] = extra_uB_proc
                        tmp_Born_PDGs.append((extra_uB_proc,isec,jsec))

            if necessary_ct_list[i*5+2] == 1:
                # Loop over sectors with final state particles only
                if isec > 2 and jsec > 2:
                    # g > g + g  
                    if id_isec == 21 and id_jsec == 21:
                        tmp_proc = 'proc_prefix_H_C_FgFg'
                    # q(qx) > g + q(qx)
                    elif id_isec == 21 and id_jsec != 21: # if there is a gluon in sector, it is always in the first position
                        tmp_proc = 'proc_prefix_H_C_FgFq'
                    # g > q(qx) + qx(q)
                    else:
                        tmp_proc = 'proc_prefix_H_C_FqFqx'

                    uB_proc = necessary_ct[i*5+2].current.shell_string_user(
                                schannel=True, forbid=True, main=False, pdg_order=False, print_id = False)
                    uB_proc_str_1 = necessary_ct[i*5+2].current.shell_string_user()
                    for j in range(0,len(uB_proc)):
                        dirpathLO = pjoin(dirpathLO_head, 'SubProcesses', "P%s" % uB_proc_str_1[j])
                        if os.path.exists(dirpathLO):
                            replace_dict_limits[tmp_proc] = uB_proc[j]
                            tmp_Born_PDGs.append((uB_proc[j],isec,jsec))
                            if uB_proc[j] not in Born_processes:
                                Born_processes.append(uB_proc[j])
                                path_Born_processes.append(dirpathLO)
                            break
                        if j == len(uB_proc) - 1:
                            extra_uB_proc = uB_proc[0]
                            replace_dict_limits[tmp_proc] = extra_uB_proc
                            tmp_Born_PDGs.append((extra_uB_proc,isec,jsec))
                    
                # Loop over sectors with at least one initial state particle
                if isec <= 2 or jsec <= 2:
                    # TODO: look at Born string for initial state singularities
                    continue

            # selection of underlying Born according to 
            # 'def compute_matrix_element_event_weight' function in ME7_integrands
            Born_PDGs.append(tmp_Born_PDGs[0])
            # write NLO_IR_limits
            filename = pjoin(dirpath, 'NLO_IR_limits_%d_%d.f' % (isec, jsec))
            file = open(pjoin(dirmadnklo,"tmp_fortran/tmp_files/NLO_IR_limits_tmp.f")).read()
            file = file % replace_dict_limits
            writer(filename).writelines(file)

        # Link LO files to each real process directory
        print('Born_processes')
        print(Born_processes)

        for i in range(0,len(Born_processes)):
            os.symlink( "%s/matrix_%s.f" % (path_Born_processes[i], Born_processes[i]), 
                        "%s/matrix_%s.f" % (dirpath, Born_processes[i]) )
            os.symlink( path_Born_processes[i] + '/%s_spin_correlations.inc' % Born_processes[i], 
                        dirpath + '/%s_spin_correlations.inc' % Born_processes[i] )
            #include all necessary leg_PDGs_* files
            for j in range(0,len(Born_PDGs)):
                if j != 0 and Born_PDGs[j][0] == Born_PDGs[j-1][0]:
                    continue
                tmp = path_Born_processes[i] + '/leg_PDGs_%s.inc' % Born_PDGs[j][0]
                if os.path.exists(tmp):
                    os.symlink( path_Born_processes[i] + '/leg_PDGs_%s.inc' % Born_PDGs[j][0], 
                        dirpath + '/leg_PDGs_%s.inc' % Born_PDGs[j][0] )


# Links to virtual dir
        for i in range(0,len(Born_processes)):
            dirpath_virtual = pjoin(dirmadnklo,glob.glob("%s/NLO_V*" % interface.user_dir_name[0])[0])
            dirpath_virtual = glob.glob("%s/SubProcesses/*%s" % (dirpath_virtual,str(Born_processes[i])))[0]
            #os.symlink(dirpath + '/virtual_recoilers.inc',dirpath_virtual + '/virtual_recoilers_%s.inc' % defining_process.shell_string(
            #                schannel=True, forbid=True, main=False, pdg_order=False, print_id = False))
            
            v_rec = recoiler_function.get_virtual_recoiler(getattr(PDGs_from_Born, "leg_PDGs_%s" % Born_processes[i]))
            data_v_rec = str(v_rec).replace('[','').replace(']','').replace(' ','').replace('(','').replace(')','')
            file = """ \
              integer, parameter :: len_iref = %d
              integer iref(2,len_iref)
              data iref/%s/ 
            """ % (len(v_rec), data_v_rec)
            filename = pjoin(dirpath_virtual, 'virtual_recoilers.inc')
            writer(filename).writelines(file)


######### Write get_Born_PDGs.f

        if len(Born_PDGs) != len(all_sector_list):
            raise MadEvent7Error('WARNING, the list of Born_PDGs is not compatible with the total number of sectors!')

        file = ''
        file += """ \
          subroutine get_Born_PDGs(isec,jsec,nexternal_Born,Born_leg_PDGs)
          implicit none
          """

        for i in range(0,len(Born_PDGs)):
            if i != 0 and Born_PDGs[i][0] == Born_PDGs[i-1][0]:
                continue
            file += """ \
                include 'leg_PDGs_%s.inc'
                """ % Born_PDGs[i][0]

        file += """ \
          integer isec, jsec
          integer nexternal_Born
          integer Born_leg_PDGs(nexternal_Born)
          \n"""


        for i in range(0,len(Born_PDGs)):
            replace_dict_tmp = {}
            replace_dict_tmp['isec'] = Born_PDGs[i][1]
            replace_dict_tmp['jsec'] = Born_PDGs[i][2]
            replace_dict_tmp['tmp_PDGs'] = 'leg_PDGS_%s' % Born_PDGs[i][0]

            if i == 0:
                replace_dict_tmp['if_elseif'] = 'if'
            else:
                replace_dict_tmp['if_elseif'] = 'elseif'

            file += """ \
               %(if_elseif)s(isec.eq.%(isec)d.and.jsec.eq.%(jsec)d) then
                  Born_leg_PDGs = %(tmp_PDGs)s \n""" % replace_dict_tmp
        
        file += """ \
          endif
          return
          end
          """

        filename = pjoin(dirpath, 'get_Born_PDGs.f')
        writer(filename).writelines(file)


######### Write driver_isec_jsec.f 

        replace_dict = {}
        for i in range(0,len(all_sector_list)):
            isec = all_sector_list[i][0]
            jsec = all_sector_list[i][1]
            replace_dict['isec'] = isec
            replace_dict['jsec'] = jsec
            # write driver
            filename = pjoin(dirpath, 'driver_%d_%d.f' % (isec, jsec))
            file = open(pjoin(dirmadnklo,"tmp_fortran/tmp_files/driver_npo_template.f")).read()
            file = file % replace_dict
            writer(filename).writelines(file)

######### Symlink for *.inc  

        os.mkdir(pjoin(dirpath, 'include'))

        path_to_math = pjoin(dirmadnklo,"Template/Fortran_tmp/src_to_common/math.inc")
        os.symlink( path_to_math, dirpath + '/include/math.inc' )

        os.symlink(dirpath + '/../../../Source/MODEL/coupl.inc',dirpath+'/include/coupl.inc')
        os.symlink(dirpath + '/../../../Source/MODEL/input.inc',dirpath+'/include/input.inc')
        os.symlink(dirpath + '/../../../Source/run.inc',dirpath+'/include/run.inc')
        os.symlink(dirpath + '/../../../Source/cuts.inc',dirpath+'/include/cuts.inc')

######### Write makefile_npo_template 

        replace_dict = {}
        proc_str = ''
        files_str = ''
        sector_str = ''
        all_str = 'all: libs'

        proc_str += """PROC_FILES= get_Born_PDGs.o matrix_%s.o""" % defining_process.shell_string(
            schannel=True, forbid=True, main=False, pdg_order=False, print_id = False)

        for i in range(0,len(Born_processes)):
            proc_str += ' matrix_' + Born_processes[i] + '.o'

        replace_dict['proc_str'] = proc_str
 
        for i in range(0,len(all_sector_list)):
            isec = all_sector_list[i][0]
            jsec = all_sector_list[i][1]
            replace_dict['isec'] = isec
            replace_dict['jsec'] = jsec
            files_str += 'FILES_%d_%d= ' % (isec, jsec)
            files_str += 'driver_%d_%d.o ' % (isec, jsec)
            files_str += 'NLO_Rsub_%d_%d.o ' % (isec, jsec)
            files_str += 'NLO_IR_limits_%d_%d.o ' % (isec, jsec)
            files_str += 'testR_%d_%d.o ' % (isec, jsec)
            files_str += 'NLO_K_%d_%d.o $(PROC_FILES) $(COMMON_FILES) \n' % (isec, jsec)
            all_str += ' sector_%d_%d' % (isec, jsec) 
            sector_str += """
sector_%d_%d_libs: libs sector_%d_%d

sector_%d_%d: $(FILES_%d_%d)
\t$(DEFAULT_F_COMPILER) $(patsubst %%,$(OBJ)/%%,$(FILES_%d_%d)) $(LIBS) -o $@ 
""" %(isec, jsec,isec, jsec,isec, jsec,isec, jsec,isec,jsec)    

        object_str = """
%.o: %.f $(INCLUDE)
\t$(DEFAULT_F_COMPILER) -c $(FFLAGS) $(FDEBUG) -o $(OBJ)/$@ $< 

%.o: $(PATH_TO_COMMON_FILES)/%.f $(INCLUDE)
\t$(DEFAULT_F_COMPILER) -c $(FFLAGS) $(FDEBUG) -o $(OBJ)/$@ $< 
"""
        replace_dict['object_str'] = object_str
        replace_dict['sector_str'] = sector_str
        replace_dict['all_str'] = all_str
        replace_dict['files_str'] = files_str
        # write makefile
        filename = pjoin(dirpath, 'makefile' )
        file = open(pjoin(dirmadnklo,"tmp_fortran/tmp_files/makefile_npo_template")).read()
        file = file % replace_dict
        writers.FileWriter(filename).write(file)



######### Write testR_template 

        replace_dict = {}
        for i in range(0,len(all_sector_list)):
            isec = all_sector_list[i][0]
            jsec = all_sector_list[i][1]
            iref = all_sector_recoilers[i]
            replace_dict['isec'] = isec
            replace_dict['jsec'] = jsec
            # replace_dict['iref'] = iref
            limit_str = ''
            is_soft = False
            is_coll = False
            list_Zsum = []
            if necessary_ct_list[i*5] == 1:
                limit_str += """
c
c     soft limit
      e1=1d0
      e2=1d0
      call do_limit_R_%d_%d(iunit,'S       ',x0,e1,e2)
"""%(isec,jsec)
                list_Zsum.append('c     call sector function ZsumSi\n')
                list_Zsum.append('        call get_Z_NLO(sNLO,sCM,alpha,%d,%d,ZsumSi,"S",ierr)\n'%(isec,jsec))
                is_soft = True
            if necessary_ct_list[i*5+1] == 1:
                limit_str += """
c
c     soft limit
      e1=1d0
      e2=1d0
      call do_limit_R_%d_%d(iunit,'S       ',x0,e1,e2)
"""%(isec,jsec)
                list_Zsum.append('c     call sector function ZsumSj\n')
                list_Zsum.append('        call get_Z_NLO(sNLO,sCM,alpha,%d,%d,ZsumSj,"S",ierr)\n'%(jsec,isec))
                is_soft = True
            # Loop over sectors with final state particles only
            if isec > 2 and jsec > 2:
                if necessary_ct_list[i*5+2] == 1:
                    limit_str += """
c
c     collinear limit
      if((iU.eq.isec.and.iS.eq.jsec.and.iB.eq.iref).or.(iU.eq.jsec.and.iS.eq.isec.and.iB.eq.iref)) then
        e1=0d0
        e2=1d0
      elseif((iU.eq.isec.and.iS.eq.iref.and.iB.eq.jsec).or.(iU.eq.jsec.and.iS.eq.iref.and.iB.eq.isec)) then
        e1=1d0
        e2=0d0
      else
        write(*,*) 'Wrong mapping for collinear limit'
        write(*,*) 'iU, iS, iB = ', iU, iS, iB
        write(*,*) 'isec, jsec = ', isec, jsec
        stop
      endif
      call do_limit_R_%d_%d(iunit,'C       ',x0,e1,e2)
"""%(isec,jsec)
                    is_coll = True
                if is_soft and is_coll:
                    limit_str += """
c
c     soft-collinear limit
      e1=e1+1d0
      e2=e2+1d0
      call do_limit_R_%d_%d(iunit,'SC      ',x0,e1,e2)
"""%(isec,jsec)
            elif isec > 2 and jsec <= 2:
                limit_str += """Collinear limits still to be specified in sectors.py """

            replace_dict['limit_str'] = limit_str
            replace_dict['NLO_proc_str'] = str(defining_process.shell_string(schannel=True, 
                                        forbid=True, main=False, pdg_order=False, print_id = False) + '_')
            str_Zsum = " ".join(list_Zsum)
            replace_dict['str_Zsum'] = str_Zsum  
            replace_dict['mapping_str'] = mapping_str
            # write testR             
            filename = pjoin(dirpath, 'testR_%d_%d.f' %(isec,jsec) )
            file = open(pjoin(dirmadnklo,"tmp_fortran/tmp_files/testR_template.f")).read()
            file = file % replace_dict
            writer(filename).writelines(file)



        return all_sectors
