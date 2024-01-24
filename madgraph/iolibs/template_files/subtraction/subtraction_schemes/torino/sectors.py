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
#import madgraph.iolibs.template_files.subtraction.subtraction_schemes.torino.sectors_for_RR as sectors_for_RR
#GIOVANNI
import stat

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
        
        # Point to the proper sector generation routine: NLO R / NNLO RR / NNLO RV
        str_contrib = contrib_definition.get_shell_name().split('_')
        if str_contrib[0] == 'NLO':
            pass
        elif str_contrib[0] == 'NNLO' and str_contrib[1] == 'RR':
            import madgraph.iolibs.template_files.subtraction.subtraction_schemes.torino.sectors_for_RR as sectors_for_RR
            all_sectors = sectors_for_RR.SectorGeneratorRR().write_RR_templates(
                                contrib_definition, defining_process, counterterms, integrated_counterterms)
            return
        elif str_contrib[0] == 'NNLO' and str_contrib[1] == 'RV':
            return #TODO RV case

        model = defining_process.get('model')
        initial_state_PDGs, final_state_PDGs = defining_process.get_cached_initial_final_pdgs()
        all_PDGs = initial_state_PDGs, final_state_PDGs

        leglist = defining_process.get('legs')
        #gl
        #print('Print from sectors.py')
        #print(defining_process.shell_string())
        #print(contrib_definition.get_shell_name())
        #print(contrib_definition.process_definition.get('id'))

        all_sectors = []
        all_sector_legs = []
        all_sector_id_legs = []
        all_sector_recoilers = []

        pert_dict = fks_common.find_pert_particles_interactions(model)
        colorlist = [model['particle_dict'][l['id']]['color'] for l in leglist]

        #print('leglist : ' + str(leglist))
        #print('colorlist : ' + str(colorlist))

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
                
                # print('list of ij : ' + str(ijlist))

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
                    #print('new_leglist NLO : ' + str(new_leglist))
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

        #print('all sector list : ' + str(all_sector_list))

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
                        if s['sector'].id[0] == 21 and s['sector'].id[1] == 21:
                            
                            #if singular_structure.substructures[0].legs[0].n == s['sector'].leg_numbers[1]:
                            #and all_legs[0].n == s['sector'].leg_numbers[1]: # should match to "j"
                            
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
        dirpath = pjoin(dirmadnklo,glob.glob("%s/NLO_R_x_R_*" % interface.user_dir_name[0])[0])
        dirpath = pjoin(dirpath, 'SubProcesses', \
                       "P%s" % defining_process.shell_string())
#        if len(glob.glob(dirpath)) == 0:
#            dirpath = ''
#            dirpath = pjoin(dirmadnklo,glob.glob("%s/NNLO_RR_x_RR_*" % interface.user_dir_name[0])[0])
#            dirpath = pjoin(dirpath, 'SubProcesses', \
#                       "P%s" % defining_process.shell_string())
            #return
#            if glob.glob(dirpath):
                #import madgraph.iolibs.template_files.subtraction.subtraction_schemes.torino.sectors_for_RR as sectors_for_RR
                #all_sectors = sectors_for_RR.SectorGeneratorRR().write_RR_templates(
                #    contrib_definition, defining_process, counterterms, integrated_counterterms,
                #            all_sectors, all_sector_legs, all_sector_id_legs, all_sector_recoilers,
                #            all_sector_list, all_sector_mass_list, all_sector_id_list,
                #            all_local_counterterms_list, necessary_ct_list, necessary_ct,
                #            dirmadnklo, dirpath)
                                                    #model, initial_state_PDGs, final_state_PDGs, all_PDGs, leglist,
                                                    #all_sectors, all_sector_legs, all_sector_id_legs, all_sector_recoilers,
                                                    #all_sector_list, all_sector_mass_list, all_sector_id_list,
                                                    #all_local_counterterms_list, necessary_ct_list, necessary_ct,
                                                    #dirmadnklo, dirpath, defining_process
                                                    #)
#                return 
#            else:
#                return
                #dirpath = ''
                #dirpath = pjoin(dirmadnklo,glob.glob("%s/NNLO_RV_x_R_*" % interface.user_dir_name[0])[0])
                #dirpath = pjoin(dirpath, 'SubProcesses', \
                #       "P%s" % defining_process.shell_string())
                #if glob.glob(dirpath):
                #    return



######### Import Born-level PDGs from proc/SupProcesses directory

        sys.path.append(pjoin(dirmadnklo,"%s/SubProcesses" % interface.user_dir_name[0]))
        import Born_PDGs as PDGs_from_Born


#==================================================================================
#  Necessary template files for real subprocess directories
# 
#       - all_sector_list.inc
#       - NLO_K_template.f
#       - NLO_Rsub_template.f
#       - driver_template.f
#       - testR_template.f
#       - NLO_IR_limits_template.f
#       - get_Born_PDGs.f
#       - makefile_npo_template.f
#       - virtual_recoiler.inc (needed for checking recoiler consistency)
#
#       - links from Born to Real subproc directories
#       - links from Born to Virtual subproc directories
# 
#==================================================================================


######### Write all_sector_list.inc

        self.write_all_sector_list_include(writers.FortranWriter, dirpath, all_sector_list)

        
######### Write NLO_K_isec_jsec.f, NLO_Rsub_isec_jsec.f

        overall_sector_info = []
        # Set replace_dict for NLO_K_isec_jsec.f
        replace_dict_ct = {}
        # Set replace_dict for NLO_Rsub_isec_jsec.f
        replace_dict_int_real = {}

        replace_dict_limits = {}
        # List of necessary underlying Born strings and particle PDGs
        Born_processes = []
        # List of dirpathLO of the necessary underlying Born
        path_Born_processes = []
        # Link LO files to each real process directory
        dirpathLO_head = pjoin(dirmadnklo,glob.glob("%s/LO_*" % interface.user_dir_name[0])[0])
        

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
            # Extract the reference particle leg from recoiler_function.py
            iref = all_sector_recoilers[i]
            replace_dict_ct['iref'] = iref
            # Check isec != jsec
            if isec == jsec:
                raise MadEvent7Error('Wrong sector indices %d,%d!' % (isec,jsec))
            replace_dict_ct['isec'] = isec
            replace_dict_ct['jsec'] = jsec
            replace_dict_limits['isec'] = isec
            replace_dict_limits['jsec'] = jsec
            replace_dict_limits['proc_prefix_real'] = str(defining_process.shell_string(schannel=True, 
                                        forbid=True, main=False, pdg_order=False, print_id = False))
            replace_dict_limits['proc_prefix_S'] = 'dummy'
            replace_dict_limits['proc_prefix_H_C_FgFg'] = 'dummy'
            replace_dict_limits['proc_prefix_H_C_FgFq'] = 'dummy'
            replace_dict_limits['proc_prefix_H_C_FqFqx'] = 'dummy'
            #replace_dict_int_real['isec'] = isec
            #replace_dict_int_real['jsec'] = jsec

            
            #replace_dict_int_real['iref'] = iref

            # Update sector_info dictionary
            sector_info = {
                'isec'          :   0,
                'jsec'          :   0,
                'iref'          :   0,
                'mapping'       :   [],
                'Born_str'      :   '',
                'Born_PDGs'     :   [],
                'path_to_Born'  :   '',
                'alt_Born_str'  :   '',
                'alt_Born_path' :   ''
            }
            sector_info['isec'] = isec
            sector_info['jsec'] = jsec
            sector_info['iref'] = iref

            


            if necessary_ct_list[i*5] == 1:
                if id_isec != 21:
                    raise MadEvent7Error('%d is not a gluon!' % isec)
                list_M2.append('if(default_soft)then\n')
                list_M2.append('KS=KS+M2_S(isec,xs,xp,wgt,ZSi,xj,xjB,nitR,1d0,wgt_chan,ierr)\n')
# ! KS=KS+M2_S2(ISEC,XS,XP,WGT,ZSI,XJ,XJB,x(1:3),NITR,1D0,wgt_chan,IERR)\n')
                list_M2.append('if(ierr.eq.1)goto 999\n')
                list_M2.append('else\n')
                list_M2.append('KS=KS+M2_S_ALT(ISEC,JSEC,IREF,XS,XP,XSB,XPB,WGT,ZSI,XJ,XJB,NITR,1D0,wgt_chan,IERR)\n')
                list_M2.append('if(ierr.eq.1)goto 999\n')
#                list_M2.append('KS=KS+M2_S_DIFF(ISEC,JSEC,IREF,XS,XP,XSB,XPB,WGT,ZSI,XJ,XJB,X(1:3),NITR,1D0,wgt_chan,IERR)\n')
#                list_M2.append('if(ierr.eq.1)goto 999\n')
                list_M2.append('endif\n')

                list_int_real.append('# call sector function ZSi\n')
                list_int_real.append('call get_Z_NLO(sNLO,sCM,alpha,isec,jsec,ZSi,"S",ierr)\n')
                list_int_real.append('if(ierr.eq.1)goto 999\n')
            if necessary_ct_list[i*5+1] == 1:
                if id_jsec != 21:
                    raise MadEvent7Error('%d is not a gluon!' % jsec)
                list_M2.append('if(default_soft)then\n')
                list_M2.append('KS=KS+M2_S(jsec,xs,xp,wgt,ZSj,xj,xjB,nitR,1d0,wgt_chan,ierr)\n')
# ! KS=KS+M2_S2(JSEC,XS,XP,WGT,ZSJ,XJ,XJB,x(1:3),NITR,1D0,wgt_chan,IERR)\n')
                list_M2.append('if(ierr.eq.1)goto 999\n')
                list_M2.append('else\n')
                list_M2.append('KS=KS+M2_S_ALT(JSEC,ISEC,IREF,XS,XP,XSB,XPB,WGT,ZSJ,XJ,XJB,NITR,1D0,wgt_chan,IERR)\n')
                list_M2.append('if(ierr.eq.1)goto 999\n')
#                list_M2.append('KS=KS+M2_S_DIFF(JSEC,ISEC,IREF,XS,XP,XSB,XPB,WGT,ZSJ,XJ,XJB,X(1:3),NITR,1D0,wgt_chan,IERR)\n')
#                list_M2.append('if(ierr.eq.1)goto 999\n')
                list_M2.append('endif\n')
                list_int_real.append('# call sector function ZSj\n')
                list_int_real.append('call get_Z_NLO(sNLO,sCM,alpha,jsec,isec,ZSj,"S",ierr)\n')
                list_int_real.append('if(ierr.eq.1)goto 999\n')
            if necessary_ct_list[i*5+2] == 1:
                # Loop over sectors with final state particles only
                if isec > 2 and jsec > 2:
                    # Check irec validity
                    if (isec == iref) or (jsec == iref):
                        raise MadEvent7Error('Wrong recoiler %d,%d,%d!' % (isec,jsec,iref))
                    # Write an identified M2_H_C_F*F* for each (**) flavour couple 
                    if id_isec == 21 and id_jsec == 21:
                        list_M2.append('KHC=KHC+M2_H_C_FgFg(isec,jsec,iref,xs,xp,xsb,xpb,wgt,xj,nitR,1d0,wgt_chan,ierr)\n')
                        list_str_defHC.append('DOUBLE PRECISION M2_H_C_FgFg')
                    elif id_isec == 21 and id_jsec != 21: # if there is a gluon in sector, it is always in the first position
                        list_M2.append('KHC=KHC+M2_H_C_FgFq(isec,jsec,iref,xs,xp,xsb,xpb,wgt,xj,nitR,1d0,wgt_chan,ierr)\n')
                        list_str_defHC.append('DOUBLE PRECISION M2_H_C_FgFq')
                    else:
                        list_M2.append('KHC=KHC+M2_H_C_FqFqx(isec,jsec,iref,xs,xp,xsb,xpb,wgt,xj,nitR,1d0,wgt_chan,ierr)\n')
                        list_str_defHC.append('DOUBLE PRECISION M2_H_C_FqFqx')
                    list_M2.append('if(ierr.eq.1)goto 999\n')

                    # default mapping for final-state collinear kernels (abc) == (ijr)
                    mapping = [('isec', isec), ('jsec', jsec), ('iref', iref)] 
                    sector_info['mapping'] = [mapping[0][1], mapping[1][1], mapping[2][1]]

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
                    iA = 1 ! default azimuth for NLO
                """ % (mapping[0][0], mapping[1][0], mapping[2][0])
                overall_sector_info.append(sector_info)
                

                str_defHC = " ".join(list_str_defHC)
                str_M2 = " ".join(list_M2)
                str_int_real = " ".join(list_int_real)
                replace_dict_ct['str_defHC'] = str_defHC
                replace_dict_ct['str_M2'] = str_M2
                replace_dict_int_real['str_int_real'] = str_int_real 
                replace_dict_int_real['NLO_process'] =  str(defining_process.shell_string(schannel=True, 
                                        forbid=True, main=False, pdg_order=False, print_id = False))
                replace_dict_int_real['mapping_str'] = mapping_str       
                replace_dict_int_real['NLO_proc_str'] = str(defining_process.shell_string(schannel=True, 
                                        forbid=True, main=False, pdg_order=False, print_id = False) + '_')
                
            if necessary_ct_list[i*5] == 1 or necessary_ct_list[i*5+1] == 1:
                # list of proc str permutations 'epem_ddx' for template
                uB_proc = necessary_ct[i*5].current.shell_string_user(
                            schannel=True, forbid=True, main=False, pdg_order=False, print_id = False)
                # list of proc str permutations '1_epem_ddx' for directory
                uB_proc_str_1 = necessary_ct[i*5].current.shell_string_user()
                for j in range(0,len(uB_proc)):
                    dirpathLO = pjoin(dirpathLO_head, 'SubProcesses', "P%s" % uB_proc_str_1[j])
                    
                    if os.path.exists(dirpathLO):
                        replace_dict_int_real['strUB'] = uB_proc[j]
                        replace_dict_limits['proc_prefix_S'] = uB_proc[j]
                        overall_sector_info[i]['Born_str'] = uB_proc[j]
                        overall_sector_info[i]['path_to_Born'] = dirpathLO
                        if uB_proc[j] not in Born_processes:
                            Born_processes.append(uB_proc[j])
                            path_Born_processes.append(dirpathLO)
                        break
                    # grouped subprocesses have no a specific LO directory
                    if j == len(uB_proc) - 1:
                        extra_uB_proc = uB_proc[0]
                        replace_dict_int_real['strUB'] = extra_uB_proc
                        replace_dict_limits['proc_prefix_S'] = extra_uB_proc
                        overall_sector_info[i]['Born_str'] = extra_uB_proc

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
                    flag = False
                    for j in range(0,len(uB_proc)):
                        dirpathLO = pjoin(dirpathLO_head, 'SubProcesses', "P%s" % uB_proc_str_1[j])                        
                        if os.path.exists(dirpathLO):
                            replace_dict_int_real['strUB'] = uB_proc[j]
                            replace_dict_limits[tmp_proc] = uB_proc[j]
                            overall_sector_info[i]['Born_str'] = uB_proc[j]
                            overall_sector_info[i]['path_to_Born'] = dirpathLO
                            if uB_proc[j] not in Born_processes:
                                Born_processes.append(uB_proc[j])
                                path_Born_processes.append(dirpathLO)
                            break
                        else:
                            list_proc = []
                            filepdg = pjoin(dirpathLO_head,'../SubProcesses/Born_PDGs.py')
                            f = open(filepdg,"r")
                            while(True):
                                line = f.readline()
                                if(line != ''):
                                    list_proc.append(line)
                                else:
                                    break 
                            f.close()
                            for k in range(len(list_proc)):
                                if(uB_proc[j] in list_proc[k]):
                                    extra_uB_proc = uB_proc[j]    
                                    replace_dict_int_real['strUB'] = extra_uB_proc
                                    replace_dict_limits[tmp_proc] = extra_uB_proc
                                    overall_sector_info[i]['Born_str'] = extra_uB_proc

                                    #gl
                                    tmp_extra_uB_proc = extra_uB_proc.split("_")
                                    fs_flavours = [x for x in tmp_extra_uB_proc[-1]]
                                    for m in range(len(fs_flavours)):
                                        if fs_flavours[m] == 's':
                                            fs_flavours[m] = 'd'
                                        elif fs_flavours[m] == 'c':
                                            fs_flavours[m] = 'u'
                                    fs_flavours = "".join(fs_flavours)
                                    tmp_extra_uB_proc[-1] = fs_flavours
                                    overall_sector_info[i]['alt_Born_str'] = "_".join(tmp_extra_uB_proc)
                                    overall_sector_info[i]['alt_Born_path'] = pjoin(dirpathLO_head, 'SubProcesses', "P%s" 
                                                                                    % "_".join(['1',overall_sector_info[i]['alt_Born_str']]))
                                    #print(overall_sector_info[i]['alt_Born_str'])
                                    #print(overall_sector_info[i]['alt_Born_path'])
                                    
                                    flag = True 
                                    break
                            if(flag == True):
                                break           
                                  
            ###   
            
            overall_sector_info[i]['Born_PDGs'] = getattr(PDGs_from_Born, "leg_PDGs_%s" % overall_sector_info[i]['Born_str'])
            # write NLO_IR_limits
            filename = pjoin(dirpath, 'NLO_IR_limits_%d_%d.f' % (isec, jsec))
            file = open(pjoin(dirmadnklo,"tmp_fortran/tmp_files/NLO_IR_limits_tmp.f")).read()
            file = file % replace_dict_limits
            writer(filename).writelines(file)

            ###

            replace_dict_int_real['isec'] = isec
            replace_dict_int_real['jsec'] = jsec
            replace_dict_int_real['iref'] = iref
            replace_dict_int_real['UBgraphs'] = overall_sector_info[i]['Born_str']
            filename_int_real = pjoin(dirpath, 'NLO_Rsub_%d_%d.f' % (isec, jsec))
            file_int_real = open(pjoin(dirmadnklo,"tmp_fortran/tmp_files/NLO_Rsub_template.f")).read()
            file_int_real = file_int_real % replace_dict_int_real
            writer(filename_int_real).writelines(file_int_real)
            UBgraphs = overall_sector_info[i]['Born_str']
            self.write_driver_npo_template(writer, dirpath, dirmadnklo, i , isec, jsec, UBgraphs)


            # write NLO_K
            filename = pjoin(dirpath, 'NLO_K_%d_%d.f' % (isec, jsec))
            file = open(pjoin(dirmadnklo,"tmp_fortran/tmp_files/NLO_K_template.f")).read()
            file = file % replace_dict_ct
            writer(filename).writelines(file)

            # write driver_npo_template
            #self.write_driver_npo_template(writer, dirpath, dirmadnklo, i , isec, jsec)

            # write testR
            self.write_testR_template_file(writer, dirpath, dirmadnklo, defining_process,
                                                    i, isec, jsec, necessary_ct_list, mapping_str)

        # check on overall_sector_info lenght
        if len(overall_sector_info) != len(all_sector_list):
            raise MadEvent7Error('WARNING, the list of sector-dictionary entries \
                                    is not compatible with the total number of sectors!')


######### Write NLO_IR_limits_isec_jsec.f and import underlying Born MEs and spin_correlations.inc

        
        
            

            

            
                

            # selection of underlying Born according to 'def compute_matrix_element_event_weight' function in ME7_integrands
            #print(dirpath)
            #print(overall_sector_info)
            
            
            


           

            #giovanni
            # write NLO_Rsub

            
            
            # for j in range(0,len(uB_proc)):
            #             dirpathLO = pjoin(dirpathLO_head, 'SubProcesses', "P%s" % uB_proc_str_1[j])
            #             if os.path.exists(dirpathLO):
            #                 replace_dict_int_real['strUB'] = uB_proc[j]
            #                 #overall_sector_info[i]['Born_str'] = uB_proc[j]
            #                 #overall_sector_info[i]['path_to_Born'] = dirpathLO
            #                 #if uB_proc[j] not in Born_processes:
            #                 #    Born_processes.append(uB_proc[j])
            #                 #    path_Born_processes.append(dirpathLO)
            #                 #break
            #             if j == len(uB_proc) - 1:
            #                 extra_uB_proc = uB_proc[0]
            #                 replace_dict_int_real['strUB'] = extra_uB_proc
                            #overall_sector_info[i]['Born_str'] = extra_uB_proc
                            


            
            
            
            
                
            #    if not overall_sector_info[i]['path_to_Born']:
            #        continue
            #    if i != 0 and overall_sector_info[i]['Born_str'] == overall_sector_info[i-1]['Born_str']:
            #        continue
            #    proc_strusr =  overall_sector_info[i]['Born_str']   
            #replace_dict_int_real['strUB'] = proc_strusr
            
            #giovanni


######### Check on real and virtual recoiler flavour

        leg_PDGs = []
        leg_PDGs.append(all_PDGs[0][0])
        leg_PDGs.append(all_PDGs[0][1])
        for i in range(0,len(final_state_PDGs)):
            leg_PDGs.append(all_PDGs[1][i])

        # Function for checking recoilers (apply get_collinear_mapped_labels, compare flavours)
        for i in range(0,len(overall_sector_info)):
            info = overall_sector_info[i]
            mapped_flavours, mapped_labels, parent_leg = recoiler_function.get_collinear_mapped_labels(
                        info['mapping'][0], info['mapping'][1], info['mapping'][2], 
                        info['isec'], info['jsec'], leg_PDGs, info['Born_PDGs']
                        )
            v_rec = recoiler_function.get_virtual_recoiler(getattr(PDGs_from_Born, "leg_PDGs_%s" % info['Born_str']))
            for j in range(0,len(v_rec)):
                Born_parent = v_rec[j][0]
                Born_recoiler = v_rec[j][1]
                if mapped_labels[parent_leg-1] == Born_parent and mapped_labels[info['iref']-1] == Born_recoiler:
                    if mapped_flavours[info['iref']-1] != info['Born_PDGs'][Born_recoiler-1]:
                        raise MadEvent7Error('Recoiler flavours from (n+1) mapping (irec = (%d,%d))    \
                                            and n virtual contribution (irec = (%d,%d)) do not match!'
                                            ) % (info['iref'], mapped_flavours[info['iref']-1], 
                                                    Born_recoiler, info['Born_PDGs'][Born_recoiler-1])


######### Write get_Born_PDGs.f

        self.write_get_Born_PDGs_file(writer, dirpath, overall_sector_info)


######### Write makefile_npo_template

        self.write_makefile_npo_file(writers.FileWriter, dirpath, dirmadnklo, defining_process, overall_sector_info)

######### Write ajob_isec_jsec

        self.write_ajob_npo_file(writers.FileWriter, dirpath, dirmadnklo, overall_sector_info)

######### Link Born files to each real process directory

        self.link_files_from_B_to_R_dir(dirpath, Born_processes, path_Born_processes, overall_sector_info)


######### Link Born files to each virtual process directory
        #self.link_files_from_B_to_V_dir(dirpath, Born_processes, path_Born_processes, Born_PDGs)


# Links to virtual dir

        for i in range(0,len(Born_processes)):
            dirpath_virtual = pjoin(dirmadnklo,glob.glob("%s/NLO_V*" % interface.user_dir_name[0])[0])
            dirpath_virtual = glob.glob("%s/SubProcesses/*%s" % (dirpath_virtual,str(Born_processes[i])))[0]

            if not glob.glob("%s/matrix.f" % dirpath_virtual):
                # symlink to Born ME
                os.symlink( "%s/matrix.f" % path_Born_processes[i], "%s/matrix.f" % dirpath_virtual )
                if len(glob.glob(dirpath_virtual + '/spin_correlations.inc')) == 0 :
                    os.symlink( path_Born_processes[i] + '/spin_correlations.inc', dirpath_virtual + '/spin_correlations.inc' )
            
                # writing virtual_recoilers.inc
                v_rec = recoiler_function.get_virtual_recoiler(getattr(PDGs_from_Born, "leg_PDGs_%s" % Born_processes[i]))
                data_v_rec = str(v_rec).replace('[','').replace(']','').replace(' ','').replace('(','').replace(')','')
                file = """ \
                  integer, parameter :: len_iref = %d
                  integer iref(2,len_iref)
                  data iref/%s/ 
                """ % (len(v_rec), data_v_rec)
                filename = pjoin(dirpath_virtual, 'virtual_recoilers.inc')
                writer(filename).writelines(file)

            if not glob.glob('%s/virtual_recoilers.inc' % dirpath):
                os.symlink( '%s/virtual_recoilers.inc' % dirpath_virtual, '%s/virtual_recoilers.inc' % dirpath)


        return all_sectors



    #===========================================================================
    # write all_sector_list include file
    #===========================================================================

    def write_all_sector_list_include(self, writer, dirpath, all_sector_list):

        replace_dict = {}
        replace_dict['len_sec_list'] = len(all_sector_list)
        replace_dict['all_sector_list'] = str(all_sector_list).replace('[','').replace(']','').replace(' ','').replace('(','').replace(')','')

        file = """ \
          integer, parameter :: lensectors = %(len_sec_list)d
          integer all_sector_list(2,lensectors)
          data all_sector_list/%(all_sector_list)s/""" % replace_dict

        filename = pjoin(dirpath, 'all_sector_list.inc')
        writer(filename).writelines(file)

        return True


    #===========================================================================
    # write driver_isec_jsec for real subprocess directory
    #===========================================================================

    def write_driver_npo_template(self, writer, dirpath, dirmadnklo, i , isec, jsec, UBgraphs):
        
        replace_dict = {}
        replace_dict['isec'] = isec
        replace_dict['jsec'] = jsec
        replace_dict['UBgraphs'] = UBgraphs

        # write driver
        filename = pjoin(dirpath, 'driver_%d_%d.f' % (isec, jsec))
        file = open(pjoin(dirmadnklo,"tmp_fortran/tmp_files/driver_npo_template.f")).read()
        file = file % replace_dict
        writer(filename).writelines(file)

        return True


    #===========================================================================
    # write file for testing limits, 'testR.f'
    #===========================================================================

    def write_testR_template_file(self, writer, dirpath, dirmadnklo, defining_process, 
                                        i, isec, jsec, necessary_ct_list, mapping_str):

        replace_dict = {}
        replace_dict['isec'] = isec
        replace_dict['jsec'] = jsec

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
            list_Zsum.append('#     call sector function ZSi\n')
            list_Zsum.append('        call get_Z_NLO(sNLO,sCM,alpha,%d,%d,ZSi,"S",ierr)\n'%(isec,jsec))
            is_soft = True
        if necessary_ct_list[i*5+1] == 1:
            limit_str += """
c
c     soft limit
      e1=0d0
      e2=1d0
      call do_limit_R_%d_%d(iunit,'S       ',x0,e1,e2)
"""%(isec,jsec)
            list_Zsum.append('#     call sector function ZSj\n')
            list_Zsum.append('        call get_Z_NLO(sNLO,sCM,alpha,%d,%d,ZSj,"S",ierr)\n'%(jsec,isec))
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
            raise MadEvent7Error('Collinear limits still to be specified in sectors.py. ')

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

        return True


    #===========================================================================
    # write 'get_Born_PDGs.f' to find labels/flavours of n-body kinematics
    #===========================================================================

    def write_get_Born_PDGs_file(self,writer, dirpath, overall_sector_info):

        file = ''
        file += """ \
          subroutine get_Born_PDGs(isec,jsec,nexternal_Born,Born_leg_PDGs)
          implicit none
          integer isec, jsec
          integer nexternal_Born
          integer Born_leg_PDGs(nexternal_Born)
          \n"""

        for i in range(0,len(overall_sector_info)):

            replace_dict_tmp = {}
            replace_dict_tmp['isec'] = overall_sector_info[i]['isec']
            replace_dict_tmp['jsec'] = overall_sector_info[i]['jsec']
            replace_dict_tmp['tmp_PDGs'] = overall_sector_info[i]['Born_PDGs']

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

        return True


    #===========================================================================
    # write 'makefile' for real subprocesses
    #===========================================================================

    def write_makefile_npo_file(self, writer, dirpath, dirmadnklo, defining_process, overall_sector_info):

        replace_dict = {}
        proc_str = ''
        files_str = ''
        sector_str = ''
        all_str = 'all: libs'
        proc_str += """PROC_FILES= get_Born_PDGs.o matrix_%s.o """ % defining_process.shell_string(
            schannel=True, forbid=True, main=False, pdg_order=False, print_id = False)
        
        for i in range(0,len(overall_sector_info)):
            if not overall_sector_info[i]['path_to_Born']:
                continue
            if i != 0 and overall_sector_info[i]['Born_str'] == overall_sector_info[i-1]['Born_str']:
                continue
            proc_str += ' matrix_' + overall_sector_info[i]['Born_str'] + '.o'

        replace_dict['proc_str'] = proc_str
 
        for i in range(0,len(overall_sector_info)):
            isec = overall_sector_info[i]['isec']
            jsec = overall_sector_info[i]['jsec']
            replace_dict['isec'] = isec
            replace_dict['jsec'] = jsec
            files_str += 'FILES_%d_%d= ' % (isec, jsec)
            files_str += 'driver_%d_%d.o ' % (isec, jsec)
            files_str += 'NLO_Rsub_%d_%d.o ' % (isec, jsec)
            files_str += 'NLO_IR_limits_%d_%d.o ' % (isec, jsec)
            if not glob.glob("%s/matrix_%s.f" % (dirpath, overall_sector_info[i]['Born_str'])):
                files_str += 'configs_%s.o ' % overall_sector_info[i]['Born_str']
                files_str += 'props_%s.o ' % overall_sector_info[i]['Born_str']
                files_str += 'decayBW_%s.o ' % overall_sector_info[i]['Born_str']
                files_str += 'leshouche_%s.o ' % overall_sector_info[i]['Born_str']

            files_str += 'testR_%d_%d.o ' % (isec, jsec)
            files_str += 'NLO_K_%d_%d.o $(PROC_FILES) $(COMMON_FILES) $(USR_FILES)\n' % (isec, jsec)
            all_str += ' sector_%d_%d' % (isec, jsec) 
            sector_str += """
sector_%d_%d_libs: libs sector_%d_%d

sector_%d_%d: $(FILES_%d_%d)
\t$(DEFAULT_F_COMPILER) $(patsubst %%,$(OBJ)/%%,$(FILES_%d_%d)) $(LIBS) $(LIBSC) -o $@ 
""" %(isec, jsec,isec, jsec,isec, jsec,isec, jsec,isec,jsec)    

        object_str = """
%.o: %.f $(INCLUDE)
\t$(DEFAULT_F_COMPILER) -c $(FFLAGS) $(FDEBUG) -o $(OBJ)/$@ $< 

#%.o: $(PATH_TO_COMMON_FILES)/%.f $(INCLUDE)
#\t$(DEFAULT_F_COMPILER) -c $(FFLAGS) $(FDEBUG) -o $(OBJ)/$@ $<

%.o: $(PATH_TO_USR_FILES)/%.f $(INCLUDE)
\t$(DEFAULT_F_COMPILER) -c $(FFLAGS) $(FDEBUG) -o $(OBJ)/$@ $<

%.o: $(PATH_TO_USR_FILES)/%.cc
\t$(DEFAULT_CPP_COMPILER) -c $(CFLAGS) $(CDEBUG) $< -o $(OBJ)/$@ $(INC)
"""
        replace_dict['object_str'] = object_str
        replace_dict['sector_str'] = sector_str
        replace_dict['all_str'] = all_str
        replace_dict['files_str'] = files_str

        # write makefile
        filename = pjoin(dirpath, 'makefile' )
        file = open(pjoin(dirmadnklo,"tmp_fortran/tmp_files/makefile_npo_template")).read()
        file = file % replace_dict
        writer(filename).write(file)

        return True



 #===========================================================================
    # write 'makefile' for real subprocesses
    #===========================================================================

    def write_ajob_npo_file(self, writer, dirpath, dirmadnklo, overall_sector_info):

        replace_dict = {}
        proc_str = ''
        sec_str = ''
        # sector_str = ''
        # all_str = 'all: libs'
        # proc_str += """PROC_FILES= get_Born_PDGs.o matrix_%s.o """ % defining_process.shell_string(
        #     schannel=True, forbid=True, main=False, pdg_order=False, print_id = False)
        
        
        
 
        for i in range(0,len(overall_sector_info)):
            isec = overall_sector_info[i]['isec']
            jsec = overall_sector_info[i]['jsec']
            #replace_dict['isec'] = isec
            #replace_dict['jsec'] = jsec
            sec_str += 'isec=%d\n' %isec
            sec_str += 'jsec=%d\n' %jsec
            sec_str += './sector_%d_%d\n' %(isec,jsec)
#             files_str += 'FILES_%d_%d= ' % (isec, jsec)
#             files_str += 'driver_%d_%d.o ' % (isec, jsec)
#             files_str += 'NLO_Rsub_%d_%d.o ' % (isec, jsec)
#             files_str += 'NLO_IR_limits_%d_%d.o ' % (isec, jsec)
#             if not glob.glob("%s/matrix_%s.f" % (dirpath, overall_sector_info[i]['Born_str'])):
#                 files_str += 'configs_%s.o ' % overall_sector_info[i]['Born_str']
#                 files_str += 'props_%s.o ' % overall_sector_info[i]['Born_str']
#                 files_str += 'decayBW_%s.o ' % overall_sector_info[i]['Born_str']
#                 files_str += 'leshouche_%s.o ' % overall_sector_info[i]['Born_str']

#             files_str += 'testR_%d_%d.o ' % (isec, jsec)
#             files_str += 'NLO_K_%d_%d.o $(PROC_FILES) $(COMMON_FILES) $(USR_FILES)\n' % (isec, jsec)
#             all_str += ' sector_%d_%d' % (isec, jsec) 
#             sector_str += """
# sector_%d_%d_libs: libs sector_%d_%d

# sector_%d_%d: $(FILES_%d_%d)
# \t$(DEFAULT_F_COMPILER) $(patsubst %%,$(OBJ)/%%,$(FILES_%d_%d)) $(LIBS) $(LIBSC) -o $@ 
# """ %(isec, jsec,isec, jsec,isec, jsec,isec, jsec,isec,jsec)    

#         object_str = """
# %.o: %.f $(INCLUDE)
# \t$(DEFAULT_F_COMPILER) -c $(FFLAGS) $(FDEBUG) -o $(OBJ)/$@ $< 

# #%.o: $(PATH_TO_COMMON_FILES)/%.f $(INCLUDE)
# #\t$(DEFAULT_F_COMPILER) -c $(FFLAGS) $(FDEBUG) -o $(OBJ)/$@ $<

# %.o: $(PATH_TO_USR_FILES)/%.f $(INCLUDE)
# \t$(DEFAULT_F_COMPILER) -c $(FFLAGS) $(FDEBUG) -o $(OBJ)/$@ $<

# %.o: $(PATH_TO_USR_FILES)/%.cc
# \t$(DEFAULT_CPP_COMPILER) -c $(CFLAGS) $(CDEBUG) $< -o $(OBJ)/$@ $(INC)
# """
#         replace_dict['object_str'] = object_str
#         replace_dict['sector_str'] = sector_str
#         replace_dict['all_str'] = all_str
            replace_dict['sec_str'] = sec_str

        # write makefile
#        filename = pjoin(dirpath, 'ajob_%d_%d' %(isec,jsec) )
        filename = pjoin(dirpath, 'ajob1')
        file = open(pjoin(dirmadnklo,"tmp_fortran/tmp_files/ajob_template")).read()
        file = file % replace_dict
        writer(filename).write(file)
        os.chmod(filename, os.stat(filename).st_mode | stat.S_IEXEC)
            
        return True






    #===========================================================================
    # function for linking files to Real subprocess directory
    #===========================================================================

    def link_files_from_B_to_R_dir(self, dirpath, Born_processes, path_Born_processes, overall_sector_info):

        for i in range(0,len(overall_sector_info)):
            
            # Set up links to additional files related to the 'Born_str', independent of flavour-dependent Born string
            if not glob.glob(dirpath + '/ngraphs_%s.inc' % overall_sector_info[i]['Born_str']):
                os.symlink(dirpath + '/../../../Common_Files/ngraphs_%s.inc' % overall_sector_info[i]['Born_str'],
                           dirpath + '/ngraphs_%s.inc' % overall_sector_info[i]['Born_str']) 
                os.symlink(dirpath + '/../../../Common_Files/configs_%s.f' % overall_sector_info[i]['Born_str'],
                           dirpath + '/configs_%s.f' % overall_sector_info[i]['Born_str'])
                os.symlink(dirpath + '/../../../Common_Files/props_%s.f' % overall_sector_info[i]['Born_str'],
                           dirpath + '/props_%s.f' % overall_sector_info[i]['Born_str'])
                os.symlink(dirpath + '/../../../Common_Files/decayBW_%s.f' % overall_sector_info[i]['Born_str'],
                           dirpath + '/decayBW_%s.f' % overall_sector_info[i]['Born_str'])
                os.symlink(dirpath + '/../../../Common_Files/leshouche_%s.f' % overall_sector_info[i]['Born_str'],
                           dirpath + '/leshouche_%s.f' % overall_sector_info[i]['Born_str'])
                
            # Set up link to matrix elements and their spin_correlation files related to the the flavour-dependent Born string
            if not overall_sector_info[i]['path_to_Born']:
                #gl
                if not glob.glob("%s/matrix_%s.f" % (dirpath, overall_sector_info[i]['alt_Born_str'])):
                    os.symlink( "%s/matrix_%s.f" % (overall_sector_info[i]['alt_Born_path'], overall_sector_info[i]['alt_Born_str']), 
                            "%s/matrix_%s.f" % (dirpath, overall_sector_info[i]['alt_Born_str']) )
                    os.symlink( overall_sector_info[i]['alt_Born_path'] + '/%s_spin_correlations.inc' % overall_sector_info[i]['alt_Born_str'], 
                            dirpath + '/%s_spin_correlations.inc' % overall_sector_info[i]['alt_Born_str'] )
                continue

            if not glob.glob("%s/matrix_%s.f" % (dirpath, overall_sector_info[i]['Born_str'])):
                os.symlink( "%s/matrix_%s.f" % (overall_sector_info[i]['path_to_Born'], overall_sector_info[i]['Born_str']), 
                            "%s/matrix_%s.f" % (dirpath, overall_sector_info[i]['Born_str']) )
                os.symlink( overall_sector_info[i]['path_to_Born'] + '/%s_spin_correlations.inc' % overall_sector_info[i]['Born_str'], 
                            dirpath + '/%s_spin_correlations.inc' % overall_sector_info[i]['Born_str'] )
            
 
                
                
                


