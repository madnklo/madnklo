#
# Functions to build necessary template files for RR subprocess directories
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
from itertools import permutations

import madgraph.iolibs.export_v4 as export_v4
import madgraph.iolibs.export_ME7 as export_ME7
import madgraph.interface.madgraph_interface as interface
import madgraph.iolibs.template_files.subtraction.subtraction_schemes.torino.sectors as sectors


class MadEvent7Error(Exception):
    pass


logger = logging.getLogger('madgraph')  

#==================================================================================
#  Necessary template files for real subprocess directories
# 
#      ? - all_sector_list.inc
#      ? - NLO_K_template.f
#      ? - NLO_Rsub_template.f
#      ? - driver_template.f
#      ? - testR_template.f
#      ? - NLO_IR_limits_template.f
#      ? - get_Born_PDGs.f
#      ? - makefile_npo_template.f
#      ? - virtual_recoiler.inc (needed for checking recoiler consistency)
#
#      ? - links from Born to Real subproc directories
#      ? - links from Born to Virtual subproc directories
# 
#==================================================================================


class SectorGeneratorRR(sectors.SectorGenerator):

    def write_RR_templates(self, contrib_definition, defining_process, counterterms, integrated_counterterms):

        model = defining_process.get('model')
        initial_state_PDGs, final_state_PDGs = defining_process.get_cached_initial_final_pdgs()
        all_PDGs = initial_state_PDGs, final_state_PDGs

        leglist = defining_process.get('legs')
        #gl
        print('INTO RR SECTOR')
        print(defining_process.shell_string())
        print(contrib_definition.get_shell_name())
        print(contrib_definition.process_definition.get('id'))

        all_sectors = []
        all_sector_legs = []
        all_sector_id_legs = []
        all_sector_recoilers = []

        pert_dict = fks_common.find_pert_particles_interactions(model)
        colorlist = [model['particle_dict'][l['id']]['color'] for l in leglist]

        # Generate sectors: in general ijkl case, i_s and k_s are the particles
        # that can go soft.
        # list of possible singular configurations, so
        # at NLO: from i and j to possible ij
        # at NNLO: from i,j,k to ijk
        #      3) [g g g], (3)
        #         [g g q], [g g bq],  (2)
        #         [g q bq],  (1)
        #         [bq q q'], [bq q bq'], [bq q q], [bq q bq] (0)
        #      4) [g g, g g],  (4)
        #         [g g, g q], [g g, g bq],  (3)
        #         [g g, bq q], [g q, g q], [g bq, g bq], [g q, g bq], [g q, g q'], [g bq, g bq'], [g q, g bq'],  (2)
        #         [g q, bq q], [g q, bq' q'],  (1)
        #         [bq q, bq q], [bq q, bq' q']  (0)

        # First divide according to the possibility of having 3 or 4 particle topologies

        ####################################################
        # 3-particle sectors
        fks_j_i = {}
        fks_k_ij = {}
        threep_sectors = []
        threep_sectors_id = []
        all_3p_sectors = []
        combs = []
        for i, col_i in zip(leglist, colorlist):
            if col_i == 1:
                continue
            if not i.get('state'):
                continue
            fks_j_i[i.get('number')] = [] # not strictly needed

            for j, col_j in zip(leglist, colorlist):
                if col_j == 1:
                    continue
                if j.get('number') == i.get('number') :
                    continue
                if not j['state']:
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

                if len(ijlist)==0:
                    continue 

                #print('i, j : ' + str(i['id']) + ', ' + str(j['id']))
                #print('3P ijlist : ' + str(ijlist))
                
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
                    #print('new_leglist ij: ' + str(new_leglist))
                    leglist_ij = new_leglist
                    if diagram_generation.Amplitude(new_process).get('diagrams'):
                        fks_j_i[i.get('number')].append(j.get('number'))
                        a_sector = {
                            'sector': None,
                            'counterterms': None,
                            'integrated_counterterms': None,
                            'recoiler' : None
                        }
                        a_sector['sector'] = sectors.Sector(leg_numbers=(i.get('number'), j.get('number')))
                        # TODO: define recoiler 
                        a_sector['recoiler'] = None
                        #a_sector['recoiler'] = recoiler_function.get_recoiler(defining_process,(i.get('number'),j.get('number')))
                        #all_sector_recoilers.append(a_sector['recoiler'].get('number'))
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

                    for k, col_k in zip(leglist, colorlist):
                        if k.get('number') == i.get('number') or k.get('number') == j.get('number'):
                            continue
                        if col_k == 1:
                            continue
                        if not k.get('state'):
                            continue                 
# gl
                        # if both i and j are gluons, then keep just the case in which i (number) < j (number)
                        if ij['id'] == 21 and k['id'] == 21 and k['state']:
                            if k.get('number') < ij.get('number') :
                                continue

                        # if j and i are quarks and antiquark in the final state, let j be the quark
                        #   this is needed in order to comply with the fct combine_ij inside fks_common
                        if ij['id'] == -k['id'] and k['state']:
                            if k['id'] < 0:
                                continue

                        fks_k_ij[ij.get('number')] = []

                        ijklist = fks_common.combine_ij(fks_common.to_fks_leg(ij, model),
                                                    fks_common.to_fks_leg(k, model),
                                                    model, pert_dict)

                        if len(ijklist)==0:
                            continue 

                        #print('ij, k : ' + str(ij['id']) + ', ' + str(k['id']))
                        #print('3P ijklist : ' + str(ijklist))

                        for ijk in ijklist:
                            # copy the defining process, remove i and j
                            # and replace them by ij.
                            new_process = copy.copy(defining_process)
                            # this is a temporary hack waiting that squared_orders for
                            #  defining_process are correctly passed
                            ##if set(new_process['squared_orders'].values()) == set([0,]):
                            # MZ this may not be optimal, but keep for the moment
                            new_process['squared_orders'] = {}
                            #print('1 : ' + str(leglist))
                            new_leglist = copy.copy(leglist_ij)
                            #print('2 : ' + str(new_leglist))
                            new_leglist[min([leglist_ij.index(ij), leglist_ij.index(k)])] = ijk
                            #print('3 : ' + str(new_leglist))
                            new_leglist.pop(max([leglist_ij.index(ij), leglist_ij.index(k)]))
                            #print('4 : ' + str(new_leglist))
                            new_process['legs'] = new_leglist
                            #print('new_leglist from ijk: ' + str(new_leglist))
                            if diagram_generation.Amplitude(new_process).get('diagrams'):
                                fks_k_ij[ij.get('number')].append(k.get('number'))
                                a_sector = {
                                    'sector': None,
                                    'counterterms': None,
                                    'integrated_counterterms': None,
                                    'recoiler' : None
                                }
                                a_sector['sector'] = sectors.Sector(leg_numbers=(ij.get('number'), k.get('number')))
                                # TODO: define recoiler 
                                a_sector['recoiler'] = None
                                #a_sector['recoiler'] = recoiler_function.get_recoiler(defining_process,(i.get('number'),j.get('number')))
                                #all_sector_recoilers.append(a_sector['recoiler'].get('number'))
                                #print('Leg number : ' + str(a_sector['sector']))
                                #gl
                                all_sector_legs.append(ij.get('number'))
                                all_sector_legs.append(k.get('number'))
                                # keep track of the masses
                                a_sector['sector'].masses = (model.get('particle_dict')[ij.get('id')]['mass'],
                                                     model.get('particle_dict')[k.get('id')]['mass'])
                                #print('Masses : ' + str(a_sector['sector'].masses))
# gl
                                # keep track of the particles' identity
                                a_sector['sector'].id = (ij.get('id'), k.get('id'))
                                all_sector_id_legs.append(ij.get('id'))
                                all_sector_id_legs.append(k.get('id'))
                                #print('Identities : ' + str(a_sector['sector'].id))

                                all_sectors.append(a_sector)


                        # TODO: remove redundant sectors; in the symmetrised case
                        # Zijkl = ijkl + ijlk + jikl + jilk + klij + klji + lkij + lkji
                        tmp_sector = [i.get('number'),j.get('number'),k.get('number')]
                        tmp_sector_id = [i.get('id'),j.get('id'),k.get('id')]

                        if len(threep_sectors) == 0 or tmp_sector not in combs:
                            a_3p_sector = {
                            'sector': None,
                            'counterterms': None,
                            'integrated_counterterms': None,
                            'recoiler' : None
                            }
                            a_3p_sector['sector'] = sectors.Sector(leg_numbers=(i.get('number'),j.get('number'),k.get('number')))
                            #a_sector['recoiler'] = recoiler_function.get_recoiler(defining_process,(i.get('number'),j.get('number')))
                            #all_3p_sector_recoilers.append(a_sector['recoiler'].get('number'))
                            # keep track of the masses
                            #a_sector['sector'].masses = (model.get('particle_dict')[i.get('id')]['mass'],
                            #                         model.get('particle_dict')[j.get('id')]['mass'])
# gl
                            # keep track of the particles' identity
                            a_3p_sector['sector'].id = (i.get('id'), j.get('id'),k.get('id'))
                            all_3p_sectors.append(a_3p_sector)
                            logger.info('NNLO sector found, legs %d, %d, %d' % a_3p_sector['sector'].leg_numbers)


                            threep_sectors.append(tmp_sector)
                            threep_sectors_id.append(tmp_sector_id)
                            combs.append([i.get('number'),j.get('number'),k.get('number')])     #ijk
                            combs.append([i.get('number'),k.get('number'),j.get('number')])     #ikj
                            combs.append([j.get('number'),i.get('number'),k.get('number')])     #jik
                            combs.append([j.get('number'),k.get('number'),i.get('number')])     #jki
                            combs.append([k.get('number'),i.get('number'),j.get('number')])     #kij
                            combs.append([k.get('number'),j.get('number'),i.get('number')])     #kji
                                
                        elif tmp_sector in combs:
                            continue

        print('3p sectors : ' + str(threep_sectors))
        print('3p sectors id. : ' + str(threep_sectors_id))



        ####################################################

        # 4-particle sectors: NLO x NLO case
        fks_j_from_i = {}
        fks_l_from_k = {}
        fourp_sectors = []
        fourp_sectors_id = []
        all_4p_sectors = []
        combs = []
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
                
                if len(ijlist)==0:
                    continue 
                
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
                        a_sector['sector'] = sectors.Sector(leg_numbers=(i.get('number'), j.get('number')))
                        # TODO: define recoiler 
                        a_sector['recoiler'] = None
                        #a_sector['recoiler'] = recoiler_function.get_recoiler(defining_process,(i.get('number'),j.get('number')))
                        #all_sector_recoilers.append(a_sector['recoiler'].get('number'))
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
                        #logger.info('First part of 4p NNLO sector found, legs %d, %d' % a_sector['sector'].leg_numbers)


                # Define k,l in [ijkl]
                for k, col_k in zip(leglist, colorlist):
                    if k.get('number') == i.get('number') or k.get('number') == j.get('number'):
                        continue
                    if col_k == 1:
                        continue
                    if not k.get('state'):
                        continue
                    fks_l_from_k[k.get('number')] = [] # not strictly needed

                    for l, col_l in zip(leglist, colorlist):
                        if col_l == 1:
                            continue
                        if l.get('number') == k.get('number') or l.get('number') == i.get('number') or l.get('number') == j.get('number'):
                            continue
                        # if k is not a gluon, then l must not be a final state gluon
                        if k['id'] != 21 and l['id'] == 21 and l['state']:
                            continue
# gl
                        # if both k and l are gluons, then keep just the case in which k (number) < l (number)
                        if k['id'] == 21 and l['id'] == 21 and l['state']:
                            if l.get('number') < k.get('number') :
                                continue

                        # if j and i are quarks and antiquark in the final state, let j be the quark
                        #   this is needed in order to comply with the fct combine_ij inside fks_common
                        if k['id'] == -l['id'] and l['state']:
                            if l['id'] < 0:
                                continue

                        kllist = fks_common.combine_ij(fks_common.to_fks_leg(k, model),
                                               fks_common.to_fks_leg(l, model),
                                               model, pert_dict)

                        if len(kllist)==0:
                            continue

                        for kl in kllist:
                            # copy the defining process, remove i and j
                            # and replace them by ij.
                            new_process = copy.copy(defining_process)
                            # this is a temporary hack waiting that squared_orders for
                            #  defining_process are correctly passed
                            ##if set(new_process['squared_orders'].values()) == set([0,]):
                            # MZ this may not be optimal, but keep for the moment
                            new_process['squared_orders'] = {}

                            new_leglist = copy.copy(leglist)
                            new_leglist[min([leglist.index(k), leglist.index(l)])] = kl
                            new_leglist.pop(max([leglist.index(k), leglist.index(l)]))
                            new_process['legs'] = new_leglist
                            if diagram_generation.Amplitude(new_process).get('diagrams'):
                                fks_l_from_k[k.get('number')].append(l.get('number'))
                                a_sector = {
                                    'sector': None,
                                    'counterterms': None,
                                    'integrated_counterterms': None,
                                    'recoiler' : None
                                }
                                a_sector['sector'] = sectors.Sector(leg_numbers=(k.get('number'), l.get('number')))
                                # TODO: define recoiler 
                                a_sector['recoiler'] = None
                                #a_sector['recoiler'] = recoiler_function.get_recoiler(defining_process,(i.get('number'),j.get('number')))
                                #all_sector_recoilers.append(a_sector['recoiler'].get('number'))
 
                                all_sector_legs.append(k.get('number'))
                                all_sector_legs.append(l.get('number'))
                                # keep track of the masses
                                a_sector['sector'].masses = (model.get('particle_dict')[k.get('id')]['mass'],
                                                     model.get('particle_dict')[l.get('id')]['mass'])
                                # keep track of the particles' identity
                                a_sector['sector'].id = (k.get('id'), l.get('id'))
                                all_sector_id_legs.append(k.get('id'))
                                all_sector_id_legs.append(l.get('id'))
                                #print('Identities : ' + str(a_sector['sector'].id))

                                all_sectors.append(a_sector)
                                #logger.info('Second part of 4p NNLO sector found, legs %d, %d' % a_sector['sector'].leg_numbers)
 

                        # TODO: remove redundant sectors; in the symmetrised case
                        # Zijkl = ijkl + ijlk + jikl + jilk + klij + klji + lkij + lkji
                        tmp_sector = [i.get('number'),j.get('number'),k.get('number'),l.get('number')]
                        tmp_sector_id = [i.get('id'),j.get('id'),k.get('id'),l.get('id')]

                        if len(fourp_sectors) == 0 or tmp_sector not in combs:
                            a_4p_sector = {
                            'sector': None,
                            'counterterms': None,
                            'integrated_counterterms': None,
                            'recoiler' : None
                            }
                            a_4p_sector['sector'] = sectors.Sector(leg_numbers=(i.get('number'),j.get('number'),k.get('number'),l.get('number')))
                            #a_sector['recoiler'] = recoiler_function.get_recoiler(defining_process,(i.get('number'),j.get('number')))
                            #all_3p_sector_recoilers.append(a_sector['recoiler'].get('number'))
                            # keep track of the masses
                            #a_sector['sector'].masses = (model.get('particle_dict')[i.get('id')]['mass'],
                            #                         model.get('particle_dict')[j.get('id')]['mass'])
# gl
                            # keep track of the particles' identity
                            a_4p_sector['sector'].id = (i.get('id'), j.get('id'),k.get('id'),l.get('id'))
                            all_4p_sectors.append(a_4p_sector)
                            logger.info('NNLO sector found, legs %d, %d, %d, %d' % a_4p_sector['sector'].leg_numbers)

                            fourp_sectors.append(tmp_sector)
                            fourp_sectors_id.append(tmp_sector_id)
                            combs.append([i.get('number'),j.get('number'),k.get('number'),l.get('number')])     #ijkl
                            combs.append([j.get('number'),i.get('number'),k.get('number'),l.get('number')])     #jikl
                            combs.append([i.get('number'),j.get('number'),l.get('number'),k.get('number')])     #ijlk
                            combs.append([j.get('number'),i.get('number'),l.get('number'),k.get('number')])     #jilk
                            combs.append([k.get('number'),l.get('number'),i.get('number'),j.get('number')])     #klij
                            combs.append([k.get('number'),l.get('number'),j.get('number'),i.get('number')])     #klji
                            combs.append([l.get('number'),k.get('number'),i.get('number'),j.get('number')])     #lkij
                            combs.append([l.get('number'),k.get('number'),j.get('number'),i.get('number')])     #lkji
                                
                        elif tmp_sector in combs:
                            continue

                        #fourp_sectors.append([i.get('number'),j.get('number'),k.get('number'),l.get('number')])
                        #fourp_sectors_id.append([i.get('id'),j.get('id'),k.get('id'),l.get('id')])

        if not all_4p_sectors:
            logger.critical('WARNING, no 4p_sectors found for %s' % defining_process.nice_string())
   
        print('All sectors for RR NNLO : ' + str(fourp_sectors))
        print('All sectors id for RR NNLO : ' + str(fourp_sectors_id))


        # Now for each sector we need to find the corresponding counterterms

        all_3p_sector_list = [s['sector'].leg_numbers for s in all_3p_sectors]  
        all_3p_sector_id_list = [s['sector'].id for s in all_3p_sectors]
        all_4p_sector_list = [s['sector'].leg_numbers for s in all_4p_sectors]  
        all_4p_sector_id_list = [s['sector'].id for s in all_4p_sectors]


        ######## 3p #########
        all_3p_local_counterterms_list = []
        all_3p_K1_ct = []
        all_3p_K2_ct = []
        all_3p_K12_ct = []
        uB_all_3p_K1_ct = []
        uB_all_3p_K2_ct = []
        for s in all_3p_sectors:
            s['sector'].all_3p_sector_list = all_3p_sector_list
            s['sector'].all_3p_sector_id_list = all_3p_sector_id_list

            if counterterms is not None:
                s['counterterms'] = []
                necessary_3p_ct1_list = [0] * (6)
                necessary_3p_ct1 = [0] * (6)
                necessary_3p_ct2_list = [0] * (10)
                necessary_3p_ct2 = [0] * (10)
                necessary_3p_ct12_list = [0] * (21)

                print('****** NEW SECTOR ******')
                print(str(s['sector'].leg_numbers[0]) + ' ' + str(s['sector'].leg_numbers[1]) + ' ' + str(s['sector'].leg_numbers[2]))
                print(str(s['sector'].id[0]) + ' ' + str(s['sector'].id[1]) + ' ' + str(s['sector'].id[2]))

                for i_ct, ct in enumerate(counterterms):
                    current = ct.nodes[0].current
                    n_subs = current.get('singular_structure').substructures
                    singular_structure = current.get('singular_structure').substructures[0]
                    all_legs = singular_structure.get_all_legs()
                    leg_numbers = s['sector'].leg_numbers
                    ileg = leg_numbers[0]
                    jleg = leg_numbers[1]
                    kleg = leg_numbers[2]

                    # safety check
                    if (not sorted([l.n for l in current.get('singular_structure').get_all_legs()]) == sorted([ileg,jleg,kleg])) and \
                        (not sorted([l.n for l in current.get('singular_structure').get_all_legs()]) == sorted([ileg,kleg])) and \
                        (not sorted([l.n for l in current.get('singular_structure').get_all_legs()]) == sorted([ileg,jleg])) and \
                        (not sorted([l.n for l in current.get('singular_structure').get_all_legs()]) == sorted([jleg,kleg])) and \
                        (not sorted([l.n for l in current.get('singular_structure').get_all_legs()]) == sorted([ileg])) and \
                        (not sorted([l.n for l in current.get('singular_structure').get_all_legs()]) == sorted([jleg])) and \
                        (not sorted([l.n for l in current.get('singular_structure').get_all_legs()]) == sorted([kleg])) :
                        continue

                    #print('GOOD CT')
                    #print('ict + ct : ' + str(i_ct) + ' ' + str(ct))
                    #print(str(ileg) + ' ' + str(jleg) + ' ' + str(kleg))
                    #print(str(s['sector'].id[0]) + ' ' + str(s['sector'].id[1]) + ' ' + str(s['sector'].id[2]))
                    # if len(current.get('singular_structure').substructures) > 1 :
                    #     print('# of subs : ' + str(len(current.get('singular_structure').substructures)))
                    #     for i in range(0,len(current.get('singular_structure').substructures)):
                    #         #print(i)
                    #         subs = current.get('singular_structure').substructures[i-1]
                    #         print('subs ' + ' : ' + str(subs))

                    #         n = len(subs.substructures) 
                    #         for j in range(1,n) :
                    #             print('subs.sub : ' + str(subs.substructures[j-1]))
                    # else :
                    #     print('singular_structure  : ' + str(singular_structure))
                    #     n = len(singular_structure.substructures) 
                    #     if n > 0 :
                    #         for k in range(0,n) :
                    #             print('subs.sub1 : ' + str(singular_structure.substructures[k-1]))
                    #             m = len(singular_structure.substructures[k-1].substructures) 
                    #             if m > 0 :
                    #                 print('subs.sub2 : ' + str(singular_structure.substructures[k-1].substructures[0]))
                                    

                    # Identify cts for L1_ijk
                    # L1_ijk  : 6  -> [Si, Sj, Sk, HCij, HCik, HCjk]    

                    if len(n_subs) == 1 and len(all_legs) == 1:
                        # Si
                        if s['sector'].id[0] == 21 :
                            s['counterterms'].append(i_ct)
                            necessary_3p_ct1_list[0] = 'S_g' 
                            necessary_3p_ct1[0] = ct
                        # Sj
                        if s['sector'].id[1] == 21 :
                            s['counterterms'].append(i_ct)
                            necessary_3p_ct1_list[1] = 'S_g' 
                            necessary_3p_ct1[1] = ct
                        # Sk
                        if s['sector'].id[2] == 21 :
                            s['counterterms'].append(i_ct)
                            necessary_3p_ct1_list[2] = 'S_g' 
                            necessary_3p_ct1[2] = ct
                    
                    if singular_structure.name()=='C' and len(all_legs)==2:
                        if not singular_structure.substructures:
                            # Cij
                            if sorted([l.n for l in all_legs]) == (sorted([ileg,jleg])):
                                s['counterterms'].append(i_ct)
                                if s['sector'].id[0] == 21 and s['sector'].id[1] == 21:
                                    necessary_3p_ct1_list[3] = 'HC_gg' 
                                elif s['sector'].id[0] == 21 and s['sector'].id[1] != 21:
                                    necessary_3p_ct1_list[3] = 'HC_gq' 
                                else :
                                    necessary_3p_ct1_list[3] = 'HC_qqx' 
                                necessary_3p_ct1[3] = ct
                            # Cik
                            if sorted([l.n for l in all_legs]) == (sorted([ileg,kleg])):
                                s['counterterms'].append(i_ct)
                                if s['sector'].id[0] == 21 and s['sector'].id[2] == 21:
                                    necessary_3p_ct1_list[4] = 'HC_gg' 
                                elif s['sector'].id[0] == 21 and s['sector'].id[2] != 21:
                                    necessary_3p_ct1_list[4] = 'HC_gq' 
                                else :
                                    necessary_3p_ct1_list[4] = 'HC_qqx' 
                                necessary_3p_ct1[4] = ct
                            # Cjk
                            if sorted([l.n for l in all_legs]) == (sorted([jleg,kleg])):
                                s['counterterms'].append(i_ct)
                                if s['sector'].id[1] == 21 and s['sector'].id[2] == 21:
                                    necessary_3p_ct1_list[5] = 'HC_gg' 
                                elif s['sector'].id[1] == 21 and s['sector'].id[2] != 21:
                                    necessary_3p_ct1_list[5] = 'HC_gq' 
                                else :
                                    necessary_3p_ct1_list[5] = 'HC_qqx' 
                                necessary_3p_ct1[5] = ct

                    # Identify cts for L2_ijk
                    # L2_ijk  : 10 -> [Sij, Sik, Sjk,
                    #                  SHCijk, SHCjik, SHCkij, 
                    #                  HCijk, 
                    #                  CijkSHCijk, CijkSHCjik, CijkSHCkij]  
                    
                    if len(n_subs) == 1 and len(singular_structure.substructures) == 0 :

                        if singular_structure.name()=='S' and len(all_legs)==2:
                            if not singular_structure.substructures:
                                # Sij
                                if sorted([l.n for l in all_legs]) == sorted([ileg,jleg]):
                                    s['counterterms'].append(i_ct)
                                    if s['sector'].id[0] == 21 and s['sector'].id[1] == 21:
                                        necessary_3p_ct2_list[0] =  'SS_gg' 
                                    else:
                                        necessary_3p_ct2_list[0] =  'SS_qqx' 
                                    necessary_3p_ct2[0] = ct
                                # Sik
                                if sorted([l.n for l in all_legs]) == sorted([ileg,kleg]):
                                    s['counterterms'].append(i_ct)
                                    if s['sector'].id[0] == 21 and s['sector'].id[2] == 21:
                                        necessary_3p_ct2_list[1] = 'SS_gg' 
                                    else:
                                        necessary_3p_ct2_list[1] = 'SS_qqx' 
                                    necessary_3p_ct2[1] = ct
                                # Sjk
                                if sorted([l.n for l in all_legs]) == sorted([jleg,kleg]):
                                    s['counterterms'].append(i_ct)
                                    if s['sector'].id[1] == 21 and s['sector'].id[2] == 21:
                                        necessary_3p_ct2_list[2] = 'SS_gg' 
                                    else:
                                        necessary_3p_ct2_list[2] = 'SS_qqx' 
                                    necessary_3p_ct2[2] = ct

                        if singular_structure.name()=='C' and len(all_legs)==3:
                            # Cijk
                            if sorted([l.n for l in all_legs]) == sorted([ileg,jleg,kleg]):
                                s['counterterms'].append(i_ct)
                                if s['sector'].id[0] == 21 and s['sector'].id[1] == 21 and s['sector'].id[2] == 21:
                                    necessary_3p_ct2_list[6] = 'HCC_ggg' 
                                elif s['sector'].id[0] == 21 and s['sector'].id[1] == 21 and s['sector'].id[2] != 21:
                                    necessary_3p_ct2_list[6] = 'HCC_ggq' 
                                elif s['sector'].id[0] == 21 and s['sector'].id[1] != 21 and  s['sector'].id[1] == (- s['sector'].id[2]):
                                    necessary_3p_ct2_list[6] = 'HCC_gqqx' 
                                else:
                                    if abs(s['sector'].id[0]) == abs(s['sector'].id[1]) and abs(s['sector'].id[1]) == abs(s['sector'].id[2]):
                                        necessary_3p_ct2_list[6] = 'HCC_qxqq' 
                                    else:
                                        necessary_3p_ct2_list[6] = 'HCC_qxqqp' 
                                necessary_3p_ct2[6] = ct

                    if len(n_subs) == 2 : 
                        # here singular_structure = coll_sub -> C(i,j)
                        all_legs_C = current.get('singular_structure').substructures[1].get_all_legs()

                        # SCijk 
                        if s['sector'].id[0] == 21 and sorted([l.n for l in all_legs_C]) == sorted([jleg,kleg]) :
                            s['counterterms'].append(i_ct)
                            if s['sector'].id[1] == 21 and s['sector'].id[2] == 21:
                                necessary_3p_ct2_list[3] = 'SHC_ggg' 
                            elif s['sector'].id[1] == 21 and s['sector'].id[2] != 21:
                                necessary_3p_ct2_list[3] = 'SHC_ggq' 
                            else:
                                necessary_3p_ct2_list[3] = 'SHC_gqqx' 
                            necessary_3p_ct2[3] = ct
                        # SCjik
                        if s['sector'].id[1] == 21 and sorted([l.n for l in all_legs_C]) == sorted([ileg,kleg]) :
                            s['counterterms'].append(i_ct)
                            if s['sector'].id[0] == 21 and s['sector'].id[2] == 21:
                                necessary_3p_ct2_list[4] = 'SHC_ggg' 
                            elif s['sector'].id[0] == 21 and s['sector'].id[2] != 21:
                                necessary_3p_ct2_list[4] = 'SHC_ggq' 
                            else:
                                necessary_3p_ct2_list[4] = 'SHC_gqqx' 
                            necessary_3p_ct2[4] = ct
                        # SCkij
                        if s['sector'].id[2] == 21 and sorted([l.n for l in all_legs_C]) == sorted([ileg,jleg]) :
                            s['counterterms'].append(i_ct)
                            if s['sector'].id[0] == 21 and s['sector'].id[1] == 21:
                                necessary_3p_ct2_list[5] = 'SHC_ggg' 
                            elif s['sector'].id[0] == 21 and s['sector'].id[1] != 21:
                                necessary_3p_ct2_list[5] = 'SHC_ggq' 
                            else:
                                necessary_3p_ct2_list[5] = 'SHC_gqqx' 
                            necessary_3p_ct2[5] = ct

                    if len(n_subs) == 1 and len(singular_structure.substructures) == 1 and \
                        singular_structure.substructures[0].name() == 'S':
                        #coll_subsub = singular_structure.substructures[0]
                        #print(coll_subsub)
                        #print(sorted([l.n for l in coll_subsub.get_all_legs()]))

                        # CijkSCijk
                        if s['sector'].id[0] == 21 : #and sorted([l.n for l in all_legs]) == sorted([jleg,kleg]) :
                            s['counterterms'].append(i_ct)
                            if s['sector'].id[1] == 21 and s['sector'].id[2] == 21:
                                necessary_3p_ct2_list[7] = 'CCSHC_ggg' 
                            elif s['sector'].id[1] == 21 and s['sector'].id[2] != 21:
                                necessary_3p_ct2_list[7] = 'CCSHC_ggq' 
                            else:
                                necessary_3p_ct2_list[7] = 'CCSHC_gqqx' 
                            necessary_3p_ct2[7] = ct
                        # CijkSCjik
                        if s['sector'].id[1] == 21 : #and sorted([l.n for l in coll_subsub.get_all_legs()]) == sorted([ileg,kleg]) :
                            s['counterterms'].append(i_ct)
                            if s['sector'].id[0] == 21 and s['sector'].id[2] == 21:
                                necessary_3p_ct2_list[8] = 'CCSHC_ggg' 
                            elif s['sector'].id[0] == 21 and s['sector'].id[2] != 21:
                                necessary_3p_ct2_list[8] = 'CCSHC_ggq'
                            else:
                                necessary_3p_ct2_list[8] = 'CCSHC_gqqx'
                            necessary_3p_ct2[8] = ct
                        # CijkSCkij
                        if s['sector'].id[2] == 21 : # and sorted([l.n for l in coll_subsub.get_all_legs()]) == sorted([ileg,jleg]) :
                            s['counterterms'].append(i_ct)
                            if s['sector'].id[0] == 21 and s['sector'].id[1] == 21:
                                necessary_3p_ct2_list[9] = 'CCSHC_ggg'
                            elif s['sector'].id[0] == 21 and s['sector'].id[1] != 21:
                                necessary_3p_ct2_list[9] = 'CCSHC_ggq'
                            else:
                                necessary_3p_ct2_list[9] = 'CCSHC_gqqx'
                            necessary_3p_ct2[9] = ct


                    # Identify cts for L12_ijk
                    # L12_ijk : 21 -> [Si Sij, Si Sik, Si SHCijk, Si HCijk',
                    #                  Sj Sij, Sj Sjk, Sj SHCjik, Sj HCijk', 
                    #                  Sk Sik, Sk Sjk, Sk SHCkij, Sk HCijk',     12              
                    #                  HCij Sij, HCij SCkij, HCij HCCijk,
                    #                  HCik Sik, HCik SCjik, HCik HCCijk,
                    #                  HCjk Sjk, HCjk SCijk, HCjk HCCijk]  9                   

                    # Si Sij
                    necessary_3p_ct12_list[0] = (''.join(('S_',necessary_3p_ct2_list[0])) \
                                                 if (necessary_3p_ct1_list[0] != 0 and necessary_3p_ct2_list[0] != 0) else 0) 
                    # Si Sik
                    necessary_3p_ct12_list[1] = (''.join(('S_',necessary_3p_ct2_list[1])) \
                                                 if (necessary_3p_ct1_list[0] != 0 and necessary_3p_ct2_list[1] != 0) else 0) 
                    # Si SHCijk = Si SCijk(1-Sij-Sik)
                    necessary_3p_ct12_list[2] = (''.join(('S_',necessary_3p_ct2_list[3])) \
                                                 if (necessary_3p_ct1_list[0] != 0 and necessary_3p_ct2_list[3] != 0) else 0)  
                    # Si HCijk' = Si Cijk(1-Sij-Sik)
                    necessary_3p_ct12_list[3] = (''.join(('S_',necessary_3p_ct2_list[6])) \
                                                 if (necessary_3p_ct1_list[0] != 0 and necessary_3p_ct2_list[6] != 0) else 0)  

                    # Sj Sij
                    necessary_3p_ct12_list[4] = (''.join(('S_',necessary_3p_ct2_list[0])) \
                                                if (necessary_3p_ct1_list[1] != 0 and necessary_3p_ct2_list[0] != 0) else 0)  
                    # Sj Sjk
                    necessary_3p_ct12_list[5] = (''.join(('S_',necessary_3p_ct2_list[2])) \
                                                if (necessary_3p_ct1_list[1] != 0 and necessary_3p_ct2_list[2] != 0) else 0) 
                    # Sj SHCjik = Sj SCjik(1-Sij-Sjk)
                    necessary_3p_ct12_list[6] = (''.join(('S_',necessary_3p_ct2_list[4])) \
                                                if (necessary_3p_ct1_list[1] != 0 and necessary_3p_ct2_list[4] != 0) else 0) 
                    # Sj HCijk' = Sj Cijk(1-Sij-Sjk)
                    necessary_3p_ct12_list[7] = (''.join(('S_',necessary_3p_ct2_list[6])) \
                                                if (necessary_3p_ct1_list[1] != 0 and necessary_3p_ct2_list[6] != 0) else 0) 

                    # Sk Sik
                    necessary_3p_ct12_list[8] = (''.join(('S_',necessary_3p_ct2_list[1])) \
                                                if (necessary_3p_ct1_list[2] != 0 and necessary_3p_ct2_list[1] != 0) else 0) 
                    # Sk Sjk
                    necessary_3p_ct12_list[9] = (''.join(('S_',necessary_3p_ct2_list[2])) \
                                                if (necessary_3p_ct1_list[2] != 0 and necessary_3p_ct2_list[2] != 0) else 0) 
                    # Sk SHCkij = Sk SCkij(1-Sik-Sjk)
                    necessary_3p_ct12_list[10] = (''.join(('S_',necessary_3p_ct2_list[5])) \
                                                if (necessary_3p_ct1_list[2] != 0 and necessary_3p_ct2_list[5] != 0) else 0) 
                    # Sk HCijk' = Sk Cijk(1-Sik-Sjk)
                    necessary_3p_ct12_list[11] = (''.join(('S_',necessary_3p_ct2_list[6])) \
                                                if (necessary_3p_ct1_list[2] != 0 and necessary_3p_ct2_list[6] != 0) else 0) 

                    # HCij Sij
                    necessary_3p_ct12_list[12] = (''.join(('HC_',necessary_3p_ct2_list[0])) \
                                                    if (necessary_3p_ct1_list[3] != 0 and necessary_3p_ct2_list[0] != 0) else 0)
                    # HCik Sik
                    necessary_3p_ct12_list[15] = (''.join(('HC_',necessary_3p_ct2_list[1])) \
                                                    if (necessary_3p_ct1_list[4] != 0 and necessary_3p_ct2_list[1] != 0) else 0) 
                    # HCjk Sjk
                    necessary_3p_ct12_list[18] = (''.join(('HC_',necessary_3p_ct2_list[2])) \
                                                    if (necessary_3p_ct1_list[5] != 0 and necessary_3p_ct2_list[2] != 0) else 0) 

                    # HCij SCkij
                    necessary_3p_ct12_list[13] = (''.join(('HC_',necessary_3p_ct2_list[5])) \
                                                    if (necessary_3p_ct1_list[3] != 0 and necessary_3p_ct2_list[5] != 0)else 0) 
                    # HCik SCjik
                    necessary_3p_ct12_list[16] = (''.join(('HC_',necessary_3p_ct2_list[4])) \
                                                    if (necessary_3p_ct1_list[4] != 0 and necessary_3p_ct2_list[4] != 0) else 0) 
                    # HCjk SCijk
                    necessary_3p_ct12_list[19] = (''.join(('HC_',necessary_3p_ct2_list[3])) \
                                                    if (necessary_3p_ct1_list[5] != 0 and necessary_3p_ct2_list[3] != 0) else 0) 
                    # HCij Cijk
                    necessary_3p_ct12_list[14] = (''.join(('HC_',necessary_3p_ct2_list[6])) \
                                                    if (necessary_3p_ct1_list[3] != 0 and necessary_3p_ct2_list[6] != 0) else 0) 
                    # HCik Cijk
                    necessary_3p_ct12_list[17] = (''.join(('HC_',necessary_3p_ct2_list[6])) \
                                                    if (necessary_3p_ct1_list[4] != 0 and necessary_3p_ct2_list[6] != 0) else 0) 
                    # HCjk Cijk 
                    necessary_3p_ct12_list[20] = (''.join(('HC_',necessary_3p_ct2_list[6])) \
                                                    if (necessary_3p_ct1_list[5] != 0 and necessary_3p_ct2_list[6] != 0) else 0) 
                    # HCij CijkSCkij
                    #necessary_3p_ct12_list[15] = (''.join(('HC_',necessary_3p_ct2_list[9])) \
                    #                                if (necessary_3p_ct1_list[3] != 0 and necessary_3p_ct2_list[9] != 0) else 0) 
                    # HCik CijkSCjik
                    #necessary_3p_ct12_list[19] = (''.join(('HC_',necessary_3p_ct2_list[8])) \
                    #                                if (necessary_3p_ct1_list[4] != 0 and necessary_3p_ct2_list[8] != 0) else 0) 
                    # HCjk CijkSCijk
                    #necessary_3p_ct12_list[23] = (''.join(('HC_',necessary_3p_ct2_list[7])) \
                    #                                if (necessary_3p_ct1_list[5] != 0 and necessary_3p_ct2_list[7] != 0) else 0) 

                print('K1 3p sector')
                print(necessary_3p_ct1_list)
                all_3p_K1_ct.append(necessary_3p_ct1_list)
                uB_all_3p_K1_ct.append(necessary_3p_ct1)

                print('K2 3p sector')
                print(necessary_3p_ct2_list)
                all_3p_K2_ct.append(necessary_3p_ct2_list)
                uB_all_3p_K2_ct.append(necessary_3p_ct2)
                    
                print('K12 3p sector')
                print(necessary_3p_ct12_list)
                all_3p_K12_ct.append(necessary_3p_ct12_list)

        ######## 3p #########
        all_4p_local_counterterms_list = []
        all_4p_K1_ct = []
        all_4p_K2_ct = []
        all_4p_K12_ct = []
        uB_all_4p_K1_ct = []
        uB_all_4p_K2_ct = []
        for s in all_4p_sectors:
            s['sector'].all_4p_sector_list = all_4p_sector_list
            s['sector'].all_4p_sector_id_list = all_4p_sector_id_list

            if counterterms is not None:
                s['counterterms'] = []
                necessary_4p_ct1_list = [0] * (6)
                necessary_4p_ct1 = [0] * (6)
                necessary_4p_ct2_list = [0] * (9)
                necessary_4p_ct2 = [0] * (9)
                necessary_4p_ct12_list = [0] * (18)

                print('****** NEW SECTOR ******')
                #print('ict + ct : ' + str(i_ct) + ' ' + str(ct))
                print(str(s['sector'].leg_numbers[0]) + ' ' + str(s['sector'].leg_numbers[1]) + ' ' + str(s['sector'].leg_numbers[2]) \
                      + ' ' + str(s['sector'].leg_numbers[3]))
                print(str(s['sector'].id[0]) + ' ' + str(s['sector'].id[1]) + ' ' + str(s['sector'].id[2]) + ' ' + str(s['sector'].id[3]))

                for i_ct, ct in enumerate(counterterms):
                    current = ct.nodes[0].current
                    n_subs = current.get('singular_structure').substructures
                    singular_structure = current.get('singular_structure').substructures[0]
                    all_legs = singular_structure.get_all_legs()
                    leg_numbers = s['sector'].leg_numbers
                    ileg = leg_numbers[0]
                    jleg = leg_numbers[1]
                    kleg = leg_numbers[2]
                    lleg = leg_numbers[3]

                    # safety check
                    if (not sorted([l.n for l in current.get('singular_structure').get_all_legs()]) == sorted([ileg,jleg,kleg,lleg])) and \
                        (not sorted([l.n for l in current.get('singular_structure').get_all_legs()]) == sorted([ileg,kleg,lleg])) and \
                        (not sorted([l.n for l in current.get('singular_structure').get_all_legs()]) == sorted([jleg,kleg,lleg])) and \
                        (not sorted([l.n for l in current.get('singular_structure').get_all_legs()]) == sorted([ileg,jleg,kleg])) and \
                        (not sorted([l.n for l in current.get('singular_structure').get_all_legs()]) == sorted([ileg,jleg,lleg])) and \
                        (not sorted([l.n for l in current.get('singular_structure').get_all_legs()]) == sorted([kleg,lleg])) and \
                        (not sorted([l.n for l in current.get('singular_structure').get_all_legs()]) == sorted([ileg,jleg])) and \
                        (not sorted([l.n for l in current.get('singular_structure').get_all_legs()]) == sorted([ileg])) and \
                        (not sorted([l.n for l in current.get('singular_structure').get_all_legs()]) == sorted([jleg])) and \
                        (not sorted([l.n for l in current.get('singular_structure').get_all_legs()]) == sorted([kleg])) and \
                        (not sorted([l.n for l in current.get('singular_structure').get_all_legs()]) == sorted([lleg])) :
                        continue

                    # print('GOOD CT')
                    # print('ict + ct : ' + str(i_ct) + ' ' + str(ct))
                    # print(str(ileg) + ' ' + str(jleg) + ' ' + str(kleg) + ' ' + str(lleg))
                    # print(str(len(n_subs)) + ', ' + str(len(singular_structure.substructures)))

                    # L1_ijkl  : 6 -> [Si, Sj, Sk, Sl, HCij, HCkl]

                    if len(n_subs) == 1 and len(all_legs) == 1:
                        # Si
                        if s['sector'].id[0] == 21 :
                            s['counterterms'].append(i_ct)
                            necessary_4p_ct1_list[0] = 'S_g' 
                            necessary_4p_ct1[0] = ct
                        # Sj
                        if s['sector'].id[1] == 21 :
                            s['counterterms'].append(i_ct)
                            necessary_4p_ct1_list[1] = 'S_g' 
                            necessary_4p_ct1[1] = ct
                        # Sk
                        if s['sector'].id[2] == 21 :
                            s['counterterms'].append(i_ct)
                            necessary_4p_ct1_list[2] = 'S_g' 
                            necessary_4p_ct1[2] = ct
                        # Sl
                        if s['sector'].id[3] == 21 :
                            s['counterterms'].append(i_ct)
                            necessary_4p_ct1_list[3] = 'S_g' 
                            necessary_4p_ct1[3] = ct

                    if singular_structure.name()=='C' and len(all_legs)==2:
                        if not singular_structure.substructures:
                            # Cij
                            if sorted([l.n for l in all_legs]) == (sorted([ileg,jleg])):
                                s['counterterms'].append(i_ct)
                                if s['sector'].id[0] == 21 and s['sector'].id[1] == 21:
                                    necessary_4p_ct1_list[4] = 'HC_gg' 
                                elif s['sector'].id[0] == 21 and s['sector'].id[1] != 21:
                                    necessary_4p_ct1_list[4] = 'HC_gq' 
                                else :
                                    necessary_4p_ct1_list[4] = 'HC_qqx' 
                                necessary_4p_ct1[4] = ct
                            # Ckl
                            if sorted([l.n for l in all_legs]) == (sorted([kleg,lleg])):
                                s['counterterms'].append(i_ct)
                                if s['sector'].id[2] == 21 and s['sector'].id[3] == 21:
                                    necessary_4p_ct1_list[5] = 'HC_gg' 
                                elif s['sector'].id[2] == 21 and s['sector'].id[3] != 21:
                                    necessary_4p_ct1_list[5] = 'HC_gq' 
                                else :
                                    necessary_4p_ct1_list[5] = 'HC_qqx' 
                                necessary_4p_ct1[5] = ct

                    # L2_ijkl  : 9  ->  [Sik, Sil, Sjk, Sjl, 
                    #                    SHCikl, SHCjkl, SHCkij, SHClij, 
                    #                    HCCijkl]
                    #                    with HCC_ijkl=Cijkl(1+Sik+Sil+Sjk+Sjl-SCikl-SCjkl-SCkij-SClij)] 

                    if len(n_subs) == 1 and len(singular_structure.substructures) == 0 :

                        if singular_structure.name()=='S' and len(all_legs)==2:
                            # Sik
                            if sorted([l.n for l in all_legs]) == sorted([ileg,kleg]):
                                s['counterterms'].append(i_ct)
                                if s['sector'].id[0] == 21 and s['sector'].id[2] == 21:
                                    necessary_4p_ct2_list[0] =  'SS_gg' 
                                else:
                                    necessary_4p_ct2_list[0] =  'SS_qqx' 
                                necessary_4p_ct2[0] = ct
                            # Sil
                            if sorted([l.n for l in all_legs]) == sorted([ileg,lleg]):
                                s['counterterms'].append(i_ct)
                                if s['sector'].id[0] == 21 and s['sector'].id[3] == 21:
                                    necessary_4p_ct2_list[1] = 'SS_gg' 
                                else:
                                    necessary_4p_ct2_list[1] = 'SS_qqx' 
                                necessary_4p_ct2[1] = ct
                            # Sjk
                            if sorted([l.n for l in all_legs]) == sorted([jleg,kleg]):
                                s['counterterms'].append(i_ct)
                                if s['sector'].id[1] == 21 and s['sector'].id[2] == 21:
                                    necessary_4p_ct2_list[2] = 'SS_gg' 
                                else:
                                    necessary_4p_ct2_list[2] = 'SS_qqx' 
                                necessary_4p_ct2[2] = ct

                            # Sjl
                            if sorted([l.n for l in all_legs]) == sorted([jleg,lleg]):
                                s['counterterms'].append(i_ct)
                                if s['sector'].id[1] == 21 and s['sector'].id[3] == 21:
                                    necessary_4p_ct2_list[3] = 'SS_gg' 
                                else:
                                    necessary_4p_ct2_list[3] = 'SS_qqx' 
                                necessary_4p_ct2[3] = ct

                    if len(n_subs) == 2 and len(singular_structure.substructures) == 0 :

                        if n_subs[0].name()=='S' and n_subs[1].name()=='C':
                            # SHCikl
                            if s['sector'].id[0] == 21 and sorted([l.n for l in n_subs[1].get_all_legs()]) == (sorted([kleg,lleg])):
                                if s['sector'].id[2] == 21 and s['sector'].id[3] == 21:
                                    necessary_4p_ct2_list[4] = 'SHC_ggg' 
                                elif s['sector'].id[2] == 21 and s['sector'].id[3] != 21:
                                    necessary_4p_ct2_list[4] = 'SHC_ggq' 
                                else:
                                    necessary_4p_ct2_list[4] = 'SHC_gqqx' 
                                necessary_4p_ct2[4] = ct
                            # SHCjkl
                            if s['sector'].id[1] == 21 and sorted([l.n for l in n_subs[1].get_all_legs()]) == (sorted([kleg,lleg])):
                                if s['sector'].id[2] == 21 and s['sector'].id[3] == 21:
                                    necessary_4p_ct2_list[5] = 'SHC_ggg' 
                                elif s['sector'].id[2] == 21 and s['sector'].id[3] != 21:
                                    necessary_4p_ct2_list[5] = 'SHC_ggq' 
                                else:
                                    necessary_4p_ct2_list[5] = 'SHC_gqqx' 
                                necessary_4p_ct2[5] = ct
                            # SHkij
                            if s['sector'].id[2] == 21 and sorted([l.n for l in n_subs[1].get_all_legs()]) == (sorted([ileg,jleg])):
                                if s['sector'].id[0] == 21 and s['sector'].id[1] == 21:
                                    necessary_4p_ct2_list[6] = 'SHC_ggg' 
                                elif s['sector'].id[0] == 21 and s['sector'].id[1] != 21:
                                    necessary_4p_ct2_list[6] = 'SHC_ggq' 
                                else:
                                    necessary_4p_ct2_list[6] = 'SHC_gqqx' 
                                necessary_4p_ct2[6] = ct
                            # SHlij
                            if s['sector'].id[3] == 21 and sorted([l.n for l in n_subs[1].get_all_legs()]) == (sorted([ileg,jleg])):
                                if s['sector'].id[0] == 21 and s['sector'].id[1] == 21:
                                    necessary_4p_ct2_list[7] = 'SHC_ggg' 
                                elif s['sector'].id[0] == 21 and s['sector'].id[1] != 21:
                                    necessary_4p_ct2_list[7] = 'SHC_ggq' 
                                else:
                                    necessary_4p_ct2_list[7] = 'SHC_gqqx' 
                                necessary_4p_ct2[7] = ct

                        # HCCijkl
                        if n_subs[0].name()=='C' and n_subs[1].name()=='C':
                            if sorted([l.n for l in n_subs[0].get_all_legs()]) == (sorted([ileg,jleg])) and \
                                sorted([l.n for l in n_subs[1].get_all_legs()]) == (sorted([kleg,lleg])) and \
                                necessary_4p_ct1_list[4] != 0 and necessary_4p_ct1_list[5] != 0:
                                if s['sector'].id[0] == 21 and s['sector'].id[1] == 21 and \
                                    s['sector'].id[2] == 21 and s['sector'].id[3] == 21:
                                    necessary_4p_ct2_list[8] = 'HCC_gggg'
                                elif s['sector'].id[0] == 21 and s['sector'].id[1] == 21 and \
                                    s['sector'].id[2] == 21 and s['sector'].id[3] != 21:
                                    necessary_4p_ct2_list[8] = 'HCC_gggq'
                                elif s['sector'].id[0] == 21 and s['sector'].id[1] == 21 and \
                                    abs(s['sector'].id[2]) == abs(s['sector'].id[3]):
                                    necessary_4p_ct2_list[8] = 'HCC_ggqxq'
                                elif s['sector'].id[0] == 21 and s['sector'].id[1] != 21 and \
                                    s['sector'].id[2] == 21 and s['sector'].id[3] != 21:
                                    necessary_4p_ct2_list[8] = 'HCC_gqgq'
                                elif s['sector'].id[0] == 21 and s['sector'].id[1] != 21 and \
                                    abs(s['sector'].id[2]) == abs(s['sector'].id[3]):
                                    necessary_4p_ct2_list[8] = 'HCC_gqqxq'
                                elif abs(s['sector'].id[0]) == abs(s['sector'].id[1]) and \
                                    abs(s['sector'].id[2]) == abs(s['sector'].id[3]):
                                    necessary_4p_ct2_list[8] = 'HCC_qxqqxq'
                                necessary_4p_ct2[8] = ct

                    # L12_ijkl : 18  ->  [Si Sik, Si Sil, Si SHCikl, 
                    #                     Sj Sjk, Sj Sjl, Sj SHCjkl,
                    #                     Sk Sik, Sk Sjk, Sk SHCkij,
                    #                     Sl Sil, Sl Sjl, Sl SHClij,
                    #                     HCij SCkij, HCij SClij, HCij HCCijkl,
                    #                     HCkl SCikl, HCkl SCjkl, HCkl HCCijkl]
 
                    if necessary_4p_ct1_list[0] != 0:
                        # Si Sik
                        necessary_4p_ct12_list[0] = (''.join(('S_',necessary_4p_ct2_list[0])) if (necessary_4p_ct2_list[0] != 0) else 0)
                        # Si Sil
                        necessary_4p_ct12_list[1] = (''.join(('S_',necessary_4p_ct2_list[1])) if (necessary_4p_ct2_list[1] != 0) else 0)
                        # Si SHCikl
                        necessary_4p_ct12_list[2] = (''.join(('S_',necessary_4p_ct2_list[4])) if (necessary_4p_ct2_list[4] != 0) else 0)

                    if necessary_4p_ct1_list[1] != 0:
                        # Sj Sjk 
                        necessary_4p_ct12_list[3] = (''.join(('S_',necessary_4p_ct2_list[2])) if (necessary_4p_ct2_list[2] != 0) else 0)
                        # Sj Sjl 
                        necessary_4p_ct12_list[4] = (''.join(('S_',necessary_4p_ct2_list[3])) if (necessary_4p_ct2_list[3] != 0) else 0)
                        # Sj SHCjkl
                        necessary_4p_ct12_list[5] = (''.join(('S_',necessary_4p_ct2_list[5])) if (necessary_4p_ct2_list[5] != 0) else 0)
                    
                    if necessary_4p_ct1_list[2] != 0:
                        # Sk Sik 
                        necessary_4p_ct12_list[6] = (''.join(('S_',necessary_4p_ct2_list[0])) if (necessary_4p_ct2_list[0] != 0) else 0)
                        # Sk Sjk 
                        necessary_4p_ct12_list[7] = (''.join(('S_',necessary_4p_ct2_list[2])) if (necessary_4p_ct2_list[2] != 0) else 0)
                        # Sk SHCkij
                        necessary_4p_ct12_list[8] = (''.join(('S_',necessary_4p_ct2_list[6])) if (necessary_4p_ct2_list[6] != 0) else 0)
                    
                    if necessary_4p_ct1_list[3] != 0:
                        # Sl Sil 
                        necessary_4p_ct12_list[9] = (''.join(('S_',necessary_4p_ct2_list[1])) if (necessary_4p_ct2_list[1] != 0) else 0)
                        # Sl Sjl 
                        necessary_4p_ct12_list[10] = (''.join(('S_',necessary_4p_ct2_list[3])) if (necessary_4p_ct2_list[3] != 0) else 0)
                        # Sl SHClij
                        necessary_4p_ct12_list[11] = (''.join(('S_',necessary_4p_ct2_list[7])) if (necessary_4p_ct2_list[7] != 0) else 0)

                    if necessary_4p_ct1_list[4] != 0:
                        # HCij SCkij 
                        necessary_4p_ct12_list[12] = (''.join(('HC_',necessary_4p_ct2_list[6])) if (necessary_4p_ct2_list[6] != 0) else 0)
                        # HCij SClij 
                        necessary_4p_ct12_list[13] = (''.join(('HC_',necessary_4p_ct2_list[7])) if (necessary_4p_ct2_list[7] != 0) else 0)
                        # HCij HCCijkl
                        necessary_4p_ct12_list[14] = (''.join(('HC_',necessary_4p_ct2_list[8])) if (necessary_4p_ct2_list[8] != 0) else 0)

                    if necessary_4p_ct1_list[5] != 0:
                        # HCkl SCikl 
                        necessary_4p_ct12_list[15] = (''.join(('HC_',necessary_4p_ct2_list[4])) if (necessary_4p_ct2_list[4] != 0) else 0)
                        # HCkl SCjkl
                        necessary_4p_ct12_list[16] = (''.join(('HC_',necessary_4p_ct2_list[5])) if (necessary_4p_ct2_list[5] != 0) else 0) 
                        # HCkl HCCijkl
                        necessary_4p_ct12_list[17] = (''.join(('HC_',necessary_4p_ct2_list[6])) if (necessary_4p_ct2_list[6] != 0) else 0)

                print('K1 4p sector')
                print(necessary_4p_ct1_list)
                all_4p_K1_ct.append(necessary_4p_ct1_list)
                uB_all_4p_K1_ct.append(necessary_4p_ct1)

                print('K2 4p sector')
                print(necessary_4p_ct2_list)
                all_4p_K2_ct.append(necessary_4p_ct2_list)
                uB_all_4p_K2_ct.append(necessary_4p_ct2)
                    
                print('K12 4p sector')
                print(necessary_4p_ct12_list)
                all_4p_K12_ct.append(necessary_4p_ct12_list)                


######### Set writer
        writer = writers.FortranWriter

        # Point to the right process directory
        dirmadnklo=os.getcwd()
        dirpath = pjoin(dirmadnklo,glob.glob("%s/NNLO_RR_x_RR_*" % interface.user_dir_name[0])[0])
        dirpath = pjoin(dirpath, 'SubProcesses', \
                      "P%s" % defining_process.shell_string())
        sys.path.append(pjoin(dirmadnklo,"%s/SubProcesses" % interface.user_dir_name[0]))
        import Born_PDGs as PDGs_from_Born
        import Real_PDGs as PDGs_from_Real

######### Write all_sector_list.inc
        self.write_all_sector_list_include(writers.FortranWriter, dirpath, all_3p_sector_list, all_4p_sector_list)

######### Useful quantities
        overall_sector_info = []
        UBorn_procs = []
        path_UBorn_procs = []
        UReal_procs = []
        path_UReal_procs = []
        dirpathB_head = pjoin(dirmadnklo,glob.glob("%s/LO_*" % interface.user_dir_name[0])[0])
        dirpathR_head = pjoin(dirmadnklo,glob.glob("%s/NLO_R_x_R_*" % interface.user_dir_name[0])[0])
            
######### Write NNLO_K_isec_jsec_ksec.f (3-particle sector)
        
        # Set replace_dict for NNLO_K
        replace_dict_ct = {}
        replace_dict_limits = {}
        replace_dict_double_real = {}
        necessary_default_3p_ct_list = ['S_g', 'HC_gg', 'HC_gq', 'HC_qqx', \
                                        'SS_gg', 'SS_qqx', \
                                        'SHC_ggg', 'SHC_ggq', 'SHC_gqqx', \
                                        'HCC_ggg', 'HCC_ggq', 'HCC_gqqx', 'HCC_qxqq', 'HCC_qxqqp', \
                                        'CCSHC_ggg', 'CCSHC_ggq', 'CCSHC_gqqx', \
                                        'S_SS_gg',  'S_SHC_ggg', 'S_SHC_ggq', 'S_SHC_gqqx', \
                                        'S_HCC_ggg', 'S_HCC_ggq', 'S_HCC_gqqx', \
                                        'HC_SS_gg', 'HC_SS_qqx', 'HC_SHC_ggg', 'HC_SHC_ggq', 'HC_SHC_gqqx', \
                                        'HC_HCC_ggg', 'HC_HCC_ggq', 'HC_HCC_gqqx', 'HC_HCC_qxqq', 'HC_HCC_qxqqp', \
                                        'HC_CCSHC_ggg', 'HC_CCSHC_ggq', 'HC_CCSHC_gqqx']
        K1_labels = ['S_i', 'S_j', 'S_k', 'HC_ij=C_ij(1-S_i-S_j)', 'HC_ik=C_ik(1-S_i-S_k)', 'HC_jk=C_jk(1-S_j-S_k)']
        K2_labels = ['SS_ij', 'SS_ik', 'SS_jk', \
                     'SHC_ijk=SC_ijk(1-SS_ij-SS_ik)', 'SHC_jik=SC_kij(1-SS_ij-SS_jk)', 'SHC_kij=SC_kij(1-SS_ik-SS_jk)', \
                     'HCC_ijk=CC_ijk(1-SS_ij-SS_jk-SS_ik)', \
                     'CCSHC_ijk=CC_ijk SC_ijk(1-SS_ij-SS_ik)', 'CCSHC_jik=CC_jik SC_jik(1-SS_ij-SS_jk)', 'CCSHC_kij=CC_kij SC_kij(1-SS_ik-SS_jk)'] 
        K12_labels = ['S_i SS_ij', 'S_i SS_ik', 'S_i SHC_ijk=S_i SC_ijk(1-SS_ij-SS_ik)', 'S_i HCC_ijk=S_i CC_ijk(1-SS_ij-SS_ik)(1-SC_ijk)', \
                      'S_j SS_ij', 'S_j SS_jk', 'S_j SHC_jik=S_j SC_jik(1-SS_ij-SS_jk)', 'S_j HCC_jik=S_j CC_jik(1-SS_ij-SS_jk)(1-SC_jik)', \
                      'S_k SS_ik', 'S_k SS_jk', 'S_k SHC_kij=S_k SC_kij(1-SS_ik-SS_jk)', 'S_k HCC_kij=S_k CC_kij(1-SS_ik-SS_jk)(1-SC_kij)', \
                      'HC_ij SS_ij=C_ij(1-S_i-S_j) SS_ij', 'HC_ij SC_kij', 'HC_ij HCC_kij=C_ij(1-S_i-S_j) CC_ijk(1-Sij-SCkij)', \
                      'HC_ik SS_ik=C_ik(1-S_i-S_k) SS_ik', 'HC_ik SC_jik', 'HC_ik HCC_jik=C_ik(1-S_i-S_k) CC_ijk(1-Sik-SCjik)', \
                      'HC_jk SS_jk=C_jk(1-S_j-S_k) SS_jk', 'HC_jk SC_ijk', 'HC_jk HCC_ijk=C_jk(1-S_j-S_k) CC_ijk(1-Sjk-SCijk)']  
        
        K1_3p_indices = ['isec', 'jsec', 'ksec', 'isec,jsec', 'isec,ksec', 'jsec,ksec']
        K2_3p_indices = ['isec,jsec', 'isec,ksec', 'jsec,ksec', \
                         'isec,jsec,ksec', 'jsec,isec,ksec', 'ksec,isec,jsec', \
                         'isec,jsec,ksec', \
                         'isec,jsec,ksec', 'jsec,isec,ksec', 'ksec,isec,jsec']
        K12_3p_indices = ['isec,jsec', 'isec,ksec', 'isec,jsec,ksec', 'isec,jsec,ksec', \
                          'isec,jsec', 'jsec,ksec', 'jsec,isec,ksec', 'jsec,isec,ksec', \
                          'isec, ksec', 'jsec,ksec', 'ksec,isec,jsec', 'ksec,isec,jsec', \
                          'isec,jsec', 'ksec,isec,jsec', 'ksec,isec,jsec', \
                          'isec,ksec', 'jsec,isec,ksec', 'jsec,isec,ksec', \
                          'jsec,ksec', 'isec,jsec,ksec', 'isec,jsec,ksec']

        for i in range(0,len(all_3p_sector_list)):
            list_str_defK1 = []
            list_str_M2_K1 = []
            list_str_defK2 = []
            list_str_M2_K2 = []
            list_str_defK12 = []
            list_str_M2_K12 = []
            mapping = []
            isec = all_3p_sector_list[i][0]
            jsec = all_3p_sector_list[i][1]
            ksec = all_3p_sector_list[i][2]
            lsec = 0
            id_isec = all_3p_sector_id_list[i][0]
            id_jsec = all_3p_sector_id_list[i][1]
            id_ksec = all_3p_sector_id_list[i][2]
            # Extract the reference particle leg from recoiler_function.py
            iref = 1 #all_sector_recoilers[i] #TODO define recoiler
            replace_dict_ct['iref'] = iref
            if (isec == iref) or (jsec == iref) or (ksec == iref):
                        raise MadEvent7Error('Wrong recoiler %d,%d,%d, %d!' % (isec,jsec,ksec,iref))
            # Check sector indicies
            if isec == jsec or isec == ksec or jsec == ksec:
                raise MadEvent7Error('Wrong sector indices %d,%d,%d!' % (isec,jsec,ksec))

            replace_dict_ct['isec'] = isec
            replace_dict_ct['jsec'] = jsec
            replace_dict_ct['ksec'] = ksec
            replace_dict_ct['lsec'] = lsec
            replace_dict_double_real['isec'] = isec
            replace_dict_double_real['jsec'] = jsec
            replace_dict_double_real['ksec'] = ksec
            replace_dict_double_real['lsec'] = lsec
            replace_dict_double_real['iref'] = iref
            replace_dict_limits['isec'] = isec
            replace_dict_limits['jsec'] = jsec
            replace_dict_limits['ksec'] = ksec
            replace_dict_limits['lsec'] = lsec
            replace_dict_limits['proc_prefix_rr'] = str(defining_process.shell_string(schannel=True, 
                                        forbid=True, main=False, pdg_order=False, print_id = False))
            replace_dict_double_real['proc_prefix_rr'] = str(defining_process.shell_string(schannel=True, 
                                       forbid=True, main=False, pdg_order=False, print_id = False) + '_')
            replace_dict_double_real['str_UBorn'] = 'dummy'
            replace_dict_double_real['UBgraphs'] = 'dummy'

            # Initialise ct routines to 'dummy'
            # TODO: remove since useless
            for k in range(0, len(necessary_default_3p_ct_list)):
                replace_dict_limits['proc_prefix_%s' % necessary_default_3p_ct_list[k]] = 'dummy'
            
            # Update sector_info dictionary
            sector_info = {
                'isec'          :   0,
                'jsec'          :   0,
                'ksec'          :   0,
                'lsec'          :   0,
                'iref'          :   0,
                'mapping'       :   [],
                'Born_str'      :   '', 
                'Born_PDGs'     :   [],
                'path_to_Born'  :   '',
                'alt_Born_str'  :   '',
                'alt_Born_path' :   '',
                'Real_str'      :   '',
                'Real_PDGs'     :   [],
                'path_to_Real'  :   '',
                'alt_Real_str'  :   '',
                'alt_Real_path' :   ''             
            }
            sector_info['isec'] = isec
            sector_info['jsec'] = jsec
            sector_info['ksec'] = ksec
            sector_info['lsec'] = lsec
            sector_info['iref'] = iref

            # default mapping
            mapping = [('isec', isec), ('jsec', jsec), ('ksec', ksec), ('lsec', lsec), ('iref', iref)] 
            sector_info['mapping'] = [mapping[0][1], mapping[1][1], mapping[2][1], mapping[3][1], mapping[4][1]]
            # specify ((isec,jsec,iref),(jsec,ksec,iref)) mapping choice
            mapping_str = """ \
                iU1 = %s 
                iS1 = %s
                iB1 = %s
                iU2 = %s
                iS2 = %s
                iB2 = %s 
                iA1 = 1 ! default azimuth for NLO
            """ % (mapping[0][0], mapping[1][0], mapping[4][0], mapping[1][0],mapping[2][0],mapping[4][0])


            # loop on K1 cts
            ct_list = []
            for j in range(0, len(all_3p_K1_ct[i])):
                if all_3p_K1_ct[i][j] ==  0:
                    continue
                if j <= 2:
                    list_str_M2_K1.append('K%s=K%s+M2_%s(%s,xs,xp,wgt,xj,xjB,nitR,1d0,wgt_chan,ierr)\n' 
                                       % (all_3p_K1_ct[i][j].split("_")[0], all_3p_K1_ct[i][j].split("_")[0], all_3p_K1_ct[i][j], K1_3p_indices[j]))
                    list_str_M2_K1.append('if(ierr.eq.1)goto 999\n')
                else:
                    list_str_M2_K1.append('K%s=K%s+M2_%s(%s,xs,xp,xsb,xpb,wgt,xj,xjB,nitR,1d0,wgt_chan,ierr)\n' 
                                       % (all_3p_K1_ct[i][j].split("_")[0], all_3p_K1_ct[i][j].split("_")[0], all_3p_K1_ct[i][j], K1_3p_indices[j]))
                    list_str_M2_K1.append('if(ierr.eq.1)goto 999\n')
                # Extract underlying real string
                self.get_uproc_str('Real', uB_all_3p_K1_ct[i][j], all_3p_K1_ct[i][j], dirpathR_head, replace_dict_limits, 
                                       replace_dict_double_real, UReal_procs, path_UReal_procs, sector_info)

                if all_3p_K1_ct[i][j] not in ct_list:
                    ct_list.append(all_3p_K1_ct[i][j])
                    tmp_str = """ 
c       %s
        DOUBLE PRECISION M2_%s""" %(K1_labels[j],all_3p_K1_ct[i][j])
                    list_str_defK1.append(tmp_str)

            # loop on K2 cts
            ct_list = []
            for j in range(0, len(all_3p_K2_ct[i])):
                if all_3p_K2_ct[i][j] ==  0:
                    continue
                if j <= 2:
                    list_str_M2_K2.append('K%s=K%s+M2_%s(%s,xs,xp,wgt,xj,xjB,nitR,1d0,wgt_chan,ierr)\n' 
                                       % (all_3p_K2_ct[i][j].split("_")[0], all_3p_K2_ct[i][j].split("_")[0], all_3p_K2_ct[i][j], K2_3p_indices[j]))
                    list_str_M2_K2.append('if(ierr.eq.1)goto 999\n')
                else:
                    list_str_M2_K2.append('K%s=K%s+M2_%s(%s,xs,xp,wgt,xj,xjB,nitR,1d0,wgt_chan,ierr)\n' 
                                       % (all_3p_K2_ct[i][j].split("_")[0], all_3p_K2_ct[i][j].split("_")[0], all_3p_K2_ct[i][j], K2_3p_indices[j]))
                    list_str_M2_K2.append('if(ierr.eq.1)goto 999\n')
                # Extract underlying Born string
                # TODO: can I have more then one underlying born x sector?
                self.get_uproc_str('Born', uB_all_3p_K2_ct[i][j], all_3p_K2_ct[i][j], dirpathB_head, replace_dict_limits, 
                                       replace_dict_double_real, UBorn_procs, path_UBorn_procs, sector_info)
                
                if all_3p_K2_ct[i][j] not in ct_list:
                    ct_list.append(all_3p_K2_ct[i][j])
                    tmp_str = """ 
c       %s
        DOUBLE PRECISION M2_%s""" %(K2_labels[j],all_3p_K2_ct[i][j])
                    list_str_defK2.append(tmp_str)

            # loop on K12 cts
            ct_list = []
            for j in range(0, len(all_3p_K12_ct[i])):
                if all_3p_K12_ct[i][j] == 0:
                    continue
                else:
                    lim = all_3p_K12_ct[i][j].split("_")[0] + '_' + all_3p_K12_ct[i][j].split("_")[1]
                if j <= 3:
                    list_str_M2_K12.append('K%s=K%s+M2_%s(isec,%s,xs,xp,wgt,xj,xjB,nitR,1d0,wgt_chan,ierr)\n' 
                                       % (lim, lim, all_3p_K12_ct[i][j], K12_3p_indices[j]))
                    list_str_M2_K12.append('if(ierr.eq.1)goto 999\n')
                elif j > 3 and j <= 7:
                    list_str_M2_K12.append('K%s=K%s+M2_%s(jsec,%s,xs,xp,wgt,xj,xjB,nitR,1d0,wgt_chan,ierr)\n' 
                                       % (lim, lim, all_3p_K12_ct[i][j], K12_3p_indices[j]))
                    list_str_M2_K12.append('if(ierr.eq.1)goto 999\n')
                elif j > 7 and j <= 11:
                    list_str_M2_K12.append('K%s=K%s+M2_%s(ksec,%s,xs,xp,wgt,xj,xjB,nitR,1d0,wgt_chan,ierr)\n' 
                                       % (lim, lim, all_3p_K12_ct[i][j], K12_3p_indices[j]))
                    list_str_M2_K12.append('if(ierr.eq.1)goto 999\n')
                elif j > 11 and j <= 14:
                    list_str_M2_K12.append('K%s=K%s+M2_%s(isec,jsec,%s,xs,xp,wgt,xj,xjB,nitR,1d0,wgt_chan,ierr)\n' 
                                       % (lim, lim, all_3p_K12_ct[i][j], K12_3p_indices[j]))
                    list_str_M2_K12.append('if(ierr.eq.1)goto 999\n')
                elif j > 14 and j <= 17:
                    list_str_M2_K12.append('K%s=K%s+M2_%s(isec,ksec,%s,xs,xp,wgt,xj,xjB,nitR,1d0,wgt_chan,ierr)\n' 
                                       % (lim, lim, all_3p_K12_ct[i][j], K12_3p_indices[j]))
                    list_str_M2_K12.append('if(ierr.eq.1)goto 999\n')
                elif j > 17 and j <= 20:
                    list_str_M2_K12.append('K%s=K%s+M2_%s(jsec,ksec,%s,xs,xp,wgt,xj,xjB,nitR,1d0,wgt_chan,ierr)\n' 
                                       % (lim, lim, all_3p_K12_ct[i][j], K12_3p_indices[j]))
                    list_str_M2_K12.append('if(ierr.eq.1)goto 999\n')
                if all_3p_K12_ct[i][j] not in ct_list:
                    ct_list.append(all_3p_K12_ct[i][j])
                    tmp_str = """
c       %s
        DOUBLE PRECISION M2_%s""" %(K12_labels[j],all_3p_K12_ct[i][j])
                    list_str_defK12.append(tmp_str)

            # update list of sector_info
            if sector_info['Born_str']:
                sector_info['Born_PDGs'] = getattr(PDGs_from_Born, "leg_PDGs_%s" % sector_info['Born_str'])
            if sector_info['Real_str']:
                sector_info['Real_PDGs'] = getattr(PDGs_from_Real, "leg_PDGs_%s" % sector_info['Real_str'])
            overall_sector_info.append(sector_info)

            # write NLO_K
            str_defK1 = " ".join(list_str_defK1)
            replace_dict_ct['str_defK1'] = str_defK1
            str_defK2 = " ".join(list_str_defK2)
            replace_dict_ct['str_defK2'] = str_defK2
            str_defK12 = " ".join(list_str_defK12)
            replace_dict_ct['str_defK12'] = str_defK12
            str_M2_K1 = " ".join(list_str_M2_K1)
            replace_dict_ct['str_M2_K1'] = str_M2_K1
            str_M2_K2 = " ".join(list_str_M2_K2)
            replace_dict_ct['str_M2_K2'] = str_M2_K2
            str_M2_K12 = " ".join(list_str_M2_K12)
            replace_dict_ct['str_M2_K12'] = str_M2_K12

            replace_dict_double_real['mapping_str'] = mapping_str

            # write NNLO_K
            filename = pjoin(dirpath, 'NNLO_K_%d_%d_%d.f' % (isec, jsec, ksec))
            file = open(pjoin(dirmadnklo,"tmp_fortran/tmp_files/NNLO_K_template.f")).read()
            file = file % replace_dict_ct
            writer(filename).writelines(file)

            # # check on sector_info
            # print('Born_str : ' + str(overall_sector_info[i]['Born_str']))
            # print('alt_Born_str : ' + str(overall_sector_info[i]['alt_Born_str']))
            # print('Born_PDGs : ' + str(overall_sector_info[i]['Born_PDGs']))
            # print('path_to_Born : ' + str(overall_sector_info[i]['path_to_Born']))
            # print('alt_Born_path : ' + str(overall_sector_info[i]['alt_Born_path']))
            # print('Real_str : ' + str(overall_sector_info[i]['Real_str']))
            # print('Real_PDGs : ' + str(overall_sector_info[i]['Real_PDGs']))
            # print('path_to_Real : ' + str(overall_sector_info[i]['path_to_Real']))
            # print('alt_Real_str : ' + str(overall_sector_info[i]['alt_Real_str']))
            # print('alt_Real_path : ' + str(overall_sector_info[i]['alt_Real_path']))
            
            # write NNLO_RRsub
            if sector_info['Born_str']:
                replace_dict_double_real['UBgraphs'] = overall_sector_info[i]['Born_str']
            else:
                if len(glob.glob("%s/ngraphs_dummy.inc" % dirpath)) == 0:
                    os.symlink(dirmadnklo + '/Template/Fortran_tmp/src_to_common/ngraphs_dummy.inc', dirpath + '/ngraphs_dummy.inc')
                    os.symlink(dirmadnklo + '/Template/Fortran_tmp/src_to_common/dummy_multich.f', dirpath + '/dummy_multich.f')
                
            filename = []
            filename = pjoin(dirpath, 'NNLO_RRsub_%d_%d_%d.f' % (isec, jsec, ksec))
            file = open(pjoin(dirmadnklo,"tmp_fortran/tmp_files/NNLO_RRsub_template.f")).read()
            file = file % replace_dict_double_real
            writer(filename).writelines(file)

            # write driver_RR

            
            self.write_driver_npt_template(writer, dirpath, dirmadnklo, i , isec, jsec, ksec, lsec, UBgraphs=None)

            self.write_testRR_3p_template_file(writer, dirpath, dirmadnklo, defining_process, 
                                    i, isec, jsec, ksec, lsec, all_3p_K1_ct, all_3p_K2_ct,all_3p_K12_ct)
            
            # write NNLO_IR_limits
            # GB TODO : test on UB strings
            NNLO_IR_limits_tmp_path = dirmadnklo + '/tmp_fortran/tmp_files/NNLO_limits/'
            filename = pjoin(dirpath, 'NNLO_IR_limits_%d_%d_%d.f' % (isec, jsec, ksec))
            file = open(NNLO_IR_limits_tmp_path + 'test_ct.f').read()
            file = file % replace_dict_limits
            writer(filename).writelines(file)
           

######### Write NNLO_K_isec_jsec_ksec_lsec.f and NNLO_R_isec_jsec_ksec_lsec (4-particle sector)
        
        # Set replace_dict
        replace_dict_ct = {}
        replace_dict_limits = {}
        replace_dict_double_real ={}

        necessary_default_4p_ct_list = ['S_g', 'HC_gg', 'HC_gq', 'HC_qqx', \
                                        'SS_gg', 'SS_qqx', \
                                        'SHC_ggg', 'SHC_ggq', 'SHC_gqqx', \
                                        'HCC_gggg', 'HCC_gggq', 'HCC_ggqqx', 'HCC_qxqqxq', 'HCC_gqgq', 'HCC_gqqxq',\
                                        'S_SS_gg',  'S_SHC_ggg', 'S_SHC_ggq', 'S_SHC_gqqx', \
                                        'HC_SHC_ggg', 'HC_SHC_ggq', 'HC_SHC_gqqx', \
                                        'HC_HCC_gggg', 'HC_HCC_gggq', 'HC_HCC_ggqqx', 'HC_HCC_qxqqxq', 'HC_HCC_gqgq', 'HC_HCC_gqqxq']
        K1_labels = ['S_i', 'S_j', 'S_k', 'S_l', 'HC_ij=C_ij(1-S_i-S_j)', 'HC_kl=C_kl(1-S_k-S_l)']
        K1_4p_indices = ['isec', 'jsec', 'ksec', 'lsec', 'isec,jsec', 'ksec,lsec']
        K2_labels = ['SS_ik', 'SS_il', 'SS_jk', 'SS_jl', \
                     'SHC_ikl=SC_ijk(1-SS_ik-SS_il)', 'SHC_jkl=SC_ijk(1-SS_jk-SS_jl)', 'SHC_kij=SC_kij(1-SS_ik-SS_jk)', 'SHC_lij=SC_lij(1-SS_il-SS_jl)', \
                     'HCC_ijkl=CC_ijkl(1+SS_ik+SS_il+SS_jk+SS_jl-SC_ikl-SC_jkl-SC_kij-SC_lij)'] 
        K2_4p_indices = ['isec,ksec', 'isec,lsec', 'jsec,ksec', 'jsec,lsec', \
                         'isec,ksec,lsec', 'jsec,ksec,lsec', 'ksec,isec,jsec', 'lsec,isec,jsec', \
                         'isec,jsec,ksec,lsec']
        K12_labels = ['S_i SS_ik', 'S_i SS_il', 'S_i SHC_ikl=S_i SC_ikl(1-SS_ik-SS_il)', \
                      'S_j SS_jk', 'S_j SS_jl', 'S_j SHC_jkl=S_j SC_jkl(1-SS_jk-SS_jl)', \
                      'S_k SS_ik', 'S_k SS_jk', 'S_k SHC_kij=S_k SC_kij(1-SS_ik-SS_jk)', \
                      'S_l SS_il', 'S_l SS_jl', 'S_l SHC_lij=S_l SC_lij(1-SS_il-SS_jl)', \
                      'HC_ij SC_kij', 'HC_ij SC_lij', 'HC_ij HCC_ijkl=C_ij(1-S_i-S_j) CC_ijkl(1-SC_kij-SC_lij)', \
                      'HC_kl SC_ikl', 'HC_kl SC_jkl', 'HC_kl HCC_ijkl=C_kl(1-S_k-S_l) CC_ijkl(1-SC_ikl-SC_jkl)']         
        K12_4p_indices = ['isec,ksec', 'isec,lsec', 'isec,ksec,lsec', \
                          'jsec,ksec', 'jsec,lsec', 'jsec,ksec,lsec', \
                          'isec,ksec', 'jsec,ksec', 'ksec,isec,jsec', \
                          'isec,lsec', 'jsec,lsec', 'lsec,isec,jsec', \
                          'ksec,isec,jsec', 'lsec,isec,jsec', 'isec,jsec,ksec,lsec', \
                          'isec,ksec,lsec', 'jsec,ksec,lsec', 'isec,jsec,ksec,lsec']
    
        for i in range(0,len(all_4p_sector_list)):
            list_str_defK1 = []
            list_str_M2_K1 = []
            list_str_defK2 = []
            list_str_M2_K2 = []
            list_str_defK12 = []
            list_str_M2_K12 = []
            mapping = []
            isec = all_4p_sector_list[i][0]
            jsec = all_4p_sector_list[i][1]
            ksec = all_4p_sector_list[i][2]
            lsec = all_4p_sector_list[i][3]
            id_isec = all_4p_sector_id_list[i][0]
            id_jsec = all_4p_sector_id_list[i][1]
            id_ksec = all_4p_sector_id_list[i][2]
            id_lsec = all_4p_sector_id_list[i][3]
            # Extract the reference particle leg from recoiler_function.py
            iref = 1 #all_sector_recoilers[i] #TODO define recoiler
            replace_dict_ct['iref'] = iref
            if (isec == iref) or (jsec == iref) or (ksec == iref) or (lsec == iref):
                        raise MadEvent7Error('Wrong recoiler %d,%d,%d,%d!' % (isec,jsec,ksec,lsec,iref))
            # Check sector indicies
            if isec == jsec or isec == ksec or jsec == ksec \
                or isec == lsec or jsec == lsec or ksec == lsec:
                raise MadEvent7Error('Wrong sector indices %d,%d,%d,%d!' % (isec,jsec,ksec,lsec))

            replace_dict_ct['isec'] = isec
            replace_dict_ct['jsec'] = jsec
            replace_dict_ct['ksec'] = ksec
            replace_dict_ct['lsec'] = lsec
            replace_dict_double_real['isec'] = isec
            replace_dict_double_real['jsec'] = jsec
            replace_dict_double_real['ksec'] = ksec
            replace_dict_double_real['lsec'] = lsec
            replace_dict_double_real['iref'] = iref
            
            replace_dict_limits['isec'] = isec
            replace_dict_limits['jsec'] = jsec
            replace_dict_limits['ksec'] = ksec
            replace_dict_limits['lsec'] = lsec
            replace_dict_limits['proc_prefix_rr'] = str(defining_process.shell_string(schannel=True, 
                                        forbid=True, main=False, pdg_order=False, print_id = False))
            replace_dict_double_real['proc_prefix_rr'] = str(defining_process.shell_string(schannel=True, 
                                        forbid=True, main=False, pdg_order=False, print_id = False) + '_')
            replace_dict_double_real['str_UBorn'] = 'dummy'
            replace_dict_double_real['UBgraphs'] = 'dummy'

            # Initialise ct routines to 'dummy'
            # TODO: remove since useless now
            for k in range(0, len(necessary_default_4p_ct_list)):
                replace_dict_limits['proc_prefix_%s' % necessary_default_4p_ct_list[k]] = 'dummy'
            
            # Update sector_info dictionary
            sector_info = {
                'isec'          :   0,
                'jsec'          :   0,
                'ksec'          :   0,
                'lsec'          :   0,
                'iref'          :   0,
                'mapping'       :   [],
                'Born_str'      :   '',
                'Born_PDGs'     :   [],
                'path_to_Born'  :   '',
                'alt_Born_str'  :   '',
                'alt_Born_path' :   '',
                'Real_str'      :   '',
                'Real_PDGs'     :   [],
                'path_to_Real'  :   '',
                'alt_Real_str'  :   '',
                'alt_Real_path' :   ''             
            }
            sector_info['isec'] = isec
            sector_info['jsec'] = jsec
            sector_info['ksec'] = ksec
            sector_info['lsec'] = lsec
            sector_info['iref'] = iref


            # default mapping
            mapping = [('isec', isec), ('jsec', jsec), ('ksec', ksec), ('lsec', lsec), ('iref', iref)] 
            sector_info['mapping'] = [mapping[0][1], mapping[1][1], mapping[2][1], mapping[3][1], mapping[4][1]]
            #specify ((isec,jsec,iref),(ksec,lsec,iref)) mapping choice
            mapping_str = """ \
                iU1 = %s 
                iS1 = %s
                iB1 = %s
                iU2 = %s
                iS2 = %s
                iB2 = %s 
                iA1 = 1 ! default azimuth for NLO
            """ % (mapping[0][0], mapping[1][0], mapping[4][0],mapping[2][0], mapping[3][0],mapping[4][0])

            replace_dict_double_real['mapping_str'] = mapping_str 
            # loop on K1 cts
            ct_list = []
            for j in range(0, len(all_4p_K1_ct[i])):
                if all_4p_K1_ct[i][j] ==  0:
                    continue
                if j <= 2:
                    list_str_M2_K1.append('K%s=K%s+M2_%s(%s,xs,xp,wgt,xj,xjB,nitR,1d0,wgt_chan,ierr)\n' 
                                       % (all_4p_K1_ct[i][j].split("_")[0], all_4p_K1_ct[i][j].split("_")[0], all_4p_K1_ct[i][j], K1_4p_indices[j]))
                    list_str_M2_K1.append('if(ierr.eq.1)goto 999\n')
                else:
                    list_str_M2_K1.append('K%s=K%s+M2_%s(%s,xs,xp,xsb,xpb,wgt,xj,xjB,nitR,1d0,wgt_chan,ierr)\n' 
                                       % (all_4p_K1_ct[i][j].split("_")[0], all_4p_K1_ct[i][j].split("_")[0], all_4p_K1_ct[i][j], K1_4p_indices[j]))
                    list_str_M2_K1.append('if(ierr.eq.1)goto 999\n')
                # Extract underlying real string        
                self.get_uproc_str('Real', uB_all_4p_K1_ct[i][j], all_4p_K1_ct[i][j], dirpathR_head, replace_dict_limits, 
                                       replace_dict_double_real, UReal_procs, path_UReal_procs, sector_info)
                
                if all_4p_K1_ct[i][j] not in ct_list:
                    ct_list.append(all_4p_K1_ct[i][j])
                    tmp_str = """ 
c       %s
        DOUBLE PRECISION M2_%s""" %(K1_labels[j],all_4p_K1_ct[i][j])
                    list_str_defK1.append(tmp_str)

            # loop on K2 cts
            ct_list = []
            for j in range(0, len(all_4p_K2_ct[i])):
                if all_4p_K2_ct[i][j] ==  0:
                    continue
                else:
                    list_str_M2_K2.append('K%s=K%s+M2_%s(%s,xs,xp,wgt,xj,xjB,nitR,1d0,wgt_chan,ierr)\n' 
                                       % (all_4p_K2_ct[i][j].split("_")[0], all_4p_K2_ct[i][j].split("_")[0], all_4p_K2_ct[i][j], K2_4p_indices[j]))
                    list_str_M2_K2.append('if(ierr.eq.1)goto 999\n')
                # Extract underlying Born string
                self.get_uproc_str('Born', uB_all_4p_K2_ct[i][j], all_4p_K2_ct[i][j], dirpathB_head, replace_dict_limits, 
                                       replace_dict_double_real, UBorn_procs, path_UBorn_procs, sector_info)
                    
                if all_4p_K2_ct[i][j] not in ct_list:
                    ct_list.append(all_4p_K2_ct[i][j])
                    tmp_str = """ 
c       %s
        DOUBLE PRECISION M2_%s""" %(K2_labels[j],all_4p_K2_ct[i][j])
                    list_str_defK2.append(tmp_str)

            # loop on K12 cts
            ct_list = []
            for j in range(0, len(all_4p_K12_ct[i])):
                if all_4p_K12_ct[i][j] == 0:
                    continue
                else:
                    lim = all_4p_K12_ct[i][j].split("_")[0] + '_' + all_4p_K12_ct[i][j].split("_")[1]
                if j <= 2:
                    list_str_M2_K12.append('K%s=K%s+M2_%s(isec,%s,xs,xp,wgt,xj,xjB,nitR,1d0,wgt_chan,ierr)\n' 
                                       % (lim, lim, all_4p_K12_ct[i][j], K12_4p_indices[j]))
                    list_str_M2_K12.append('if(ierr.eq.1)goto 999\n')
                elif j > 2 and j <= 5:
                    list_str_M2_K12.append('K%s=K%s+M2_%s(jsec,%s,xs,xp,wgt,xj,xjB,nitR,1d0,wgt_chan,ierr)\n' 
                                       % (lim, lim, all_4p_K12_ct[i][j], K12_4p_indices[j]))
                    list_str_M2_K12.append('if(ierr.eq.1)goto 999\n')
                elif j > 5 and j <= 8:
                    list_str_M2_K12.append('K%s=K%s+M2_%s(ksec,%s,xs,xp,wgt,xj,xjB,nitR,1d0,wgt_chan,ierr)\n' 
                                       % (lim, lim, all_4p_K12_ct[i][j], K12_4p_indices[j]))
                    list_str_M2_K12.append('if(ierr.eq.1)goto 999\n')
                elif j > 8 and j <= 11:
                    list_str_M2_K12.append('K%s=K%s+M2_%s(lsec,%s,xs,xp,wgt,xj,xjB,nitR,1d0,wgt_chan,ierr)\n' 
                                       % (lim, lim, all_4p_K12_ct[i][j], K12_4p_indices[j]))
                    list_str_M2_K12.append('if(ierr.eq.1)goto 999\n')
                elif j > 11 and j <= 14:
                    list_str_M2_K12.append('K%s=K%s+M2_%s(isec,jsec,%s,xs,xp,wgt,xj,xjB,nitR,1d0,wgt_chan,ierr)\n' 
                                       % (lim, lim, all_4p_K12_ct[i][j], K12_4p_indices[j]))
                    list_str_M2_K12.append('if(ierr.eq.1)goto 999\n')
                elif j > 14 and j <= 17:
                    list_str_M2_K12.append('K%s=K%s+M2_%s(ksec,lsec,%s,xs,xp,wgt,xj,xjB,nitR,1d0,wgt_chan,ierr)\n' 
                                       % (lim, lim, all_4p_K12_ct[i][j], K12_4p_indices[j]))
                    list_str_M2_K12.append('if(ierr.eq.1)goto 999\n')
                if all_4p_K12_ct[i][j] not in ct_list:
                    ct_list.append(all_4p_K12_ct[i][j])
                    tmp_str = """
c       %s
        DOUBLE PRECISION M2_%s""" %(K12_labels[j],all_4p_K12_ct[i][j])
                    list_str_defK12.append(tmp_str)

            # update list of sector_info
            if sector_info['Born_str']:
                sector_info['Born_PDGs'] = getattr(PDGs_from_Born, "leg_PDGs_%s" % sector_info['Born_str'])
            if sector_info['Real_str']:
                sector_info['Real_PDGs'] = getattr(PDGs_from_Real, "leg_PDGs_%s" % sector_info['Real_str'])
            overall_sector_info.append(sector_info)

            # define dictionary for NNLO_K
            str_defK1 = " ".join(list_str_defK1)
            replace_dict_ct['str_defK1'] = str_defK1
            str_defK2 = " ".join(list_str_defK2)
            replace_dict_ct['str_defK2'] = str_defK2
            str_defK12 = " ".join(list_str_defK12)
            replace_dict_ct['str_defK12'] = str_defK12
            str_M2_K1 = " ".join(list_str_M2_K1)
            replace_dict_ct['str_M2_K1'] = str_M2_K1
            str_M2_K2 = " ".join(list_str_M2_K2)
            replace_dict_ct['str_M2_K2'] = str_M2_K2
            str_M2_K12 = " ".join(list_str_M2_K12)
            replace_dict_ct['str_M2_K12'] = str_M2_K12

            # write NNLO_K
            filename = pjoin(dirpath, 'NNLO_K_%d_%d_%d_%d.f' % (isec, jsec, ksec, lsec))
            file = open(pjoin(dirmadnklo,"tmp_fortran/tmp_files/NNLO_K_template.f")).read()
            file = file % replace_dict_ct
            writer(filename).writelines(file)

            # # check on sector_info
            # print('Born_str : ' + str(overall_sector_info[i+len(all_3p_sector_list)]['Born_str']))
            # print('alt_Born_str : ' + str(overall_sector_info[i+len(all_3p_sector_list)]['alt_Born_str']))
            # print('Born_PDGs : ' + str(overall_sector_info[i+len(all_3p_sector_list)]['Born_PDGs']))
            # print('path_to_Born : ' + str(overall_sector_info[i+len(all_3p_sector_list)]['path_to_Born']))
            # print('alt_Born_path : ' + str(overall_sector_info[i+len(all_3p_sector_list)]['alt_Born_path']))
            # print('Real_str : ' + str(overall_sector_info[i+len(all_3p_sector_list)]['Real_str']))
            # print('Real_PDGs : ' + str(overall_sector_info[i+len(all_3p_sector_list)]['Real_PDGs']))
            # print('path_to_Real : ' + str(overall_sector_info[i+len(all_3p_sector_list)]['path_to_Real']))
            # print('alt_Real_str : ' + str(overall_sector_info[i+len(all_3p_sector_list)]['alt_Real_str']))
            # print('alt_Real_path : ' + str(overall_sector_info[i+len(all_3p_sector_list)]['alt_Real_path']))
            # print('  ')

            # write NNLO_RRsub
            if sector_info['Born_str']:
                replace_dict_double_real['UBgraphs'] = overall_sector_info[i+len(all_3p_sector_list)]['Born_str']
            else:
                if len(glob.glob("%s/ngraphs_dummy.inc" % dirpath)) == 0:
                    os.symlink(dirmadnklo + '/Template/Fortran_tmp/src_to_common/ngraphs_dummy.inc', dirpath + '/ngraphs_dummy.inc')
                    os.symlink(dirmadnklo + '/Template/Fortran_tmp/src_to_common/dummy_multich.f', dirpath + '/dummy_multich.f')

            filename = []
            filename = pjoin(dirpath, 'NNLO_RRsub_%d_%d_%d_%d.f' % (isec, jsec, ksec,lsec))
            file = open(pjoin(dirmadnklo,"tmp_fortran/tmp_files/NNLO_RRsub_template.f")).read()
            file = file % replace_dict_double_real
            writer(filename).writelines(file)
           
            self.write_driver_npt_template(writer, dirpath, dirmadnklo, i , isec, jsec, ksec, lsec, UBgraphs=None)


            self.write_testRR_4p_template_file(writer, dirpath, dirmadnklo, defining_process, 
                                    i, isec, jsec, ksec, lsec,all_4p_K1_ct, all_4p_K2_ct,all_4p_K12_ct)
           

#---------- Functions outside loop on sectors ----------#
######### Write get_Born_PDGs.f & get_Real_PDGs.f

        self.write_get_Born_PDGs_file(writer, dirpath, overall_sector_info)
        self.write_get_Real_PDGs_file(writer, dirpath, overall_sector_info)
        self.write_get_UnderLying_PDGs_file(writer, dirpath, overall_sector_info)

######### Write makefile_npt_template
        

######### Write ajob_isec_jsec_ksec_lsec
        

######### Link Born & Real files to each real process directory

        self.link_files_to_RR_dir(dirpath, overall_sector_info)

            
        return #all_sectors






    #===========================================================================
    # get the underlying Real/Born process strings
    #===========================================================================

    def get_uproc_str(self, u_str, ct, ct_name, dirpath, replace_dict_limits, replace_dict_double_real, 
                      proc_dir, path_proc_dir, sector_info):

        UProc = ct.current.shell_string_user(
                schannel=True, forbid=True, main=False, pdg_order=False, print_id = False)
        UProc_1 = ct.current.shell_string_user()

        flag = False
        for i in range(0,len(UProc)):
            dirpathUProc = pjoin(dirpath, 'SubProcesses', "P%s" % UProc_1[i])                    
            if os.path.exists(dirpathUProc):
                replace_dict_double_real['str_U%s' % u_str] = UProc[i]
                replace_dict_limits['proc_prefix_%s' % ct_name] = UProc[i]
                sector_info['%s_str' % u_str] = UProc[i]
                sector_info['path_to_%s' % u_str] = dirpathUProc
                if UProc[i] not in proc_dir:
                    proc_dir.append(UProc[i])
                    path_proc_dir.append(dirpathUProc)
                    break
            # if subproc has no associated directory
            else:
                list_proc = []
                filepdg = pjoin(dirpath,'../SubProcesses/%s_PDGs.py' % u_str)
                f = open(filepdg,"r")
                while(True):
                    line = f.readline()
                    if(line != ''):
                        list_proc.append(line)
                    else:
                        break 
                f.close()
                for k in range(len(list_proc)):
                    if(UProc[i] in list_proc[k]):
                        extra_UProc = UProc[i]    
                        replace_dict_double_real['str_U%s' % u_str] = UProc[i]
                        replace_dict_limits['proc_prefix_%s' % ct_name] = UProc[i]
                        sector_info['%s_str' % u_str] = UProc[i]

                        tmp_extra_UProc = extra_UProc.split("_")
                        fs_flavours = [x for x in tmp_extra_UProc[-1]]
                        for m in range(len(fs_flavours)):
                            if fs_flavours[m] == 's':
                                fs_flavours[m] = 'd'
                            elif fs_flavours[m] == 'c':
                                fs_flavours[m] = 'u'
                        fs_flavours = "".join(fs_flavours)
                        tmp_extra_UProc[-1] = fs_flavours
                        sector_info['alt_%s_str' % u_str] = "_".join(tmp_extra_UProc)
                        sector_info['alt_%s_path' % u_str] = pjoin(dirpath, 'SubProcesses', "P%s" 
                                                                                    % "_".join(['1', sector_info['alt_%s_str' % u_str]]))
                                    
                        flag = True 
                        break
                if(flag == True):
                    break 

            # if no specific underlying directory
            #if i == len(UProc) - 1:
            #    replace_dict_double_real['str_U%s' % u_str] = UProc[0]
            #    replace_dict_limits['proc_prefix_%s' % ct_name] = UProc[0]
            #    overall_sector_info['%s_str' % u_str] = UProc[0]


    #===========================================================================
    # write all_sector_list include file
    #===========================================================================

    def write_all_sector_list_include(self, writer, dirpath, all_3p_sector_list, all_4p_sector_list):

        # Define 3p sectors as 4p sectors with lsec = 0
        all_3p_sector_list_with_0 = []
        for i in range(0,len(all_3p_sector_list)):
            new_list = list(all_3p_sector_list[i])
            new_list.append(0)
            new_list = tuple(new_list)
            all_3p_sector_list_with_0.append(new_list)
        # Check
        if len(all_3p_sector_list_with_0) != len(all_3p_sector_list):
            print('WARNING: Wrong number of 3p sectors!')
            return
        
        all_sector_list = all_3p_sector_list_with_0 + all_4p_sector_list

        replace_dict = {}
        replace_dict['len_sec_list'] = len(all_sector_list)
        replace_dict['all_sector_list'] = str(all_sector_list).replace('[','').replace(']','').replace(' ','').replace('(','').replace(')','')

        file = """ \
          integer, parameter :: lensectors = %(len_sec_list)d
          integer all_sector_list(4,lensectors)
          data all_sector_list/%(all_sector_list)s/""" % replace_dict
        
        filename = pjoin(dirpath, 'all_sector_list.inc')
        writer(filename).writelines(file)

        return True
    


    #===========================================================================
    # write 'get_Born_PDGs.f' & 'get_Real_PDGs.f' to find labels/flavours of n(+1)-body kinematics
    #===========================================================================

    def write_get_Born_PDGs_file(self, writer, dirpath, overall_sector_info):

        file = ''
        file += """ \
          subroutine get_Born_PDGs(isec,jsec,ksec,lsec,nexternal_Born,Born_leg_PDGs)
          implicit none
          integer isec, jsec, ksec, lsec
          integer nexternal_Born
          integer Born_leg_PDGs(nexternal_Born)
          \n"""

        for i in range(0,len(overall_sector_info)):

            replace_dict_tmp = {}
            replace_dict_tmp['isec'] = overall_sector_info[i]['isec']
            replace_dict_tmp['jsec'] = overall_sector_info[i]['jsec']
            replace_dict_tmp['ksec'] = overall_sector_info[i]['ksec']
            replace_dict_tmp['lsec'] = overall_sector_info[i]['lsec']
            replace_dict_tmp['tmp_PDGs'] = overall_sector_info[i]['Born_PDGs']

            if i == 0:
                replace_dict_tmp['if_elseif'] = 'if'
            else:
                replace_dict_tmp['if_elseif'] = 'elseif'

            file += """ \
               %(if_elseif)s(isec.eq.%(isec)d.and.jsec.eq.%(jsec)d.and.ksec.eq.%(ksec)d.and.lsec.eq.%(lsec)d) then
                  Born_leg_PDGs = %(tmp_PDGs)s \n""" % replace_dict_tmp
        
        file += """ \
          endif
          return
          end
          """

        filename = pjoin(dirpath, 'get_Born_PDGs.f')
        writer(filename).writelines(file)

        return True
    
    def write_get_Real_PDGs_file(self, writer, dirpath, overall_sector_info):

        file = ''
        file += """ \
          subroutine get_Real_PDGs(isec,jsec,ksec,lsec,nexternal_Real,Real_leg_PDGs)
          implicit none
          integer isec, jsec, ksec, lsec
          integer nexternal_Real
          integer Real_leg_PDGs(nexternal_Real)
          \n"""

        for i in range(0,len(overall_sector_info)):

            replace_dict_tmp = {}
            replace_dict_tmp['isec'] = overall_sector_info[i]['isec']
            replace_dict_tmp['jsec'] = overall_sector_info[i]['jsec']
            replace_dict_tmp['ksec'] = overall_sector_info[i]['ksec']
            replace_dict_tmp['lsec'] = overall_sector_info[i]['lsec']
            replace_dict_tmp['tmp_PDGs'] = overall_sector_info[i]['Real_PDGs']

            if i == 0:
                replace_dict_tmp['if_elseif'] = 'if'
            else:
                replace_dict_tmp['if_elseif'] = 'elseif'

            file += """ \
               %(if_elseif)s(isec.eq.%(isec)d.and.jsec.eq.%(jsec)d.and.ksec.eq.%(ksec)d.and.lsec.eq.%(lsec)d) then
                  Real_leg_PDGs = %(tmp_PDGs)s \n""" % replace_dict_tmp
        
        file += """ \
          endif
          return
          end
          """

        filename = pjoin(dirpath, 'get_Real_PDGs.f')
        writer(filename).writelines(file)

        return True
    

    def write_get_UnderLying_PDGs_file(self, writer, dirpath, overall_sector_info):

        file = ''
        file += """ \
          subroutine get_UnderLying_PDGs(isec,jsec,ksec,lsec,npart,Underlying_leg_PDGs)
          implicit none
          include 'nexternal.inc'
          integer isec, jsec, ksec, lsec
          integer npart
          integer Underlying_leg_PDGs(npart)
          integer Born_leg_PDGs(nexternal_Born), Real_leg_PDGs(nexternal_Real) 
          \n"""

        for i in range(0,len(overall_sector_info)):

            replace_dict_tmp = {}
            replace_dict_tmp['isec'] = overall_sector_info[i]['isec']
            replace_dict_tmp['jsec'] = overall_sector_info[i]['jsec']
            replace_dict_tmp['ksec'] = overall_sector_info[i]['ksec']
            replace_dict_tmp['lsec'] = overall_sector_info[i]['lsec']
            replace_dict_tmp['tmp_Real_PDGs'] = overall_sector_info[i]['Real_PDGs']
            replace_dict_tmp['tmp_Born_PDGs'] = overall_sector_info[i]['Born_PDGs']

            if i == 0:
                replace_dict_tmp['if_elseif'] = 'if'
            else:
                replace_dict_tmp['if_elseif'] = 'elseif'

            file += """ \
               %(if_elseif)s(isec.eq.%(isec)d.and.jsec.eq.%(jsec)d.and.ksec.eq.%(ksec)d.and.lsec.eq.%(lsec)d) then
                  Real_leg_PDGs = %(tmp_Real_PDGs)s
                  Born_leg_PDGs = %(tmp_Born_PDGs)s \n""" % replace_dict_tmp
        
        file += """ \
          endif
          if(npart .eq. nexternal_Real) then
            Underlying_leg_PDGs = Real_leg_PDGs
          elseif(npart .eq. nexternal_Born) then
            Underlying_leg_PDGs = Born_leg_PDGs
          else
            write(*,*) 'Get_Underlying_PDGs: error'
            write(*,*) 'npart is neither equal to nexternal_Real nor to nexternal_Born:'
            write(*,*) 'npart, nexternal_Real, nexternal_Born = ',  npart, nexternal_Real, nexternal_Born
            stop
          endif        
          return
          end
          """

        filename = pjoin(dirpath, 'get_UnderlyingProc_PDGs.f')
        writer(filename).writelines(file)

        return True
    

    #===========================================================================
    # function for linking files to RR subprocess directory
    #===========================================================================

    def link_files_to_RR_dir(self, dirpath, overall_sector_info):

        for i in range(0,len(overall_sector_info)):

            # copy born stuff
            if not overall_sector_info[i]['Born_str']:
                continue
            
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
                
            # copy real stuff
            if not overall_sector_info[i]['Real_str']:
                continue
            if not overall_sector_info[i]['path_to_Real']:
                if not glob.glob("%s/matrix_%s.f" % (dirpath, overall_sector_info[i]['alt_Real_str'])):
                    os.symlink( "%s/matrix_%s.f" % (overall_sector_info[i]['alt_Real_path'], overall_sector_info[i]['alt_Real_str']), 
                            "%s/matrix_%s.f" % (dirpath, overall_sector_info[i]['alt_Real_str']) )
                    os.symlink( overall_sector_info[i]['alt_Real_path'] + '/%s_spin_correlations.inc' % overall_sector_info[i]['alt_Real_str'], 
                            dirpath + '/%s_spin_correlations.inc' % overall_sector_info[i]['alt_Real_str'] )
                continue

            if not glob.glob("%s/matrix_%s.f" % (dirpath, overall_sector_info[i]['Real_str'])):
                os.symlink( "%s/matrix_%s.f" % (overall_sector_info[i]['path_to_Real'], overall_sector_info[i]['Real_str']), 
                            "%s/matrix_%s.f" % (dirpath, overall_sector_info[i]['Real_str']) )
                os.symlink( overall_sector_info[i]['path_to_Real'] + '/%s_spin_correlations.inc' % overall_sector_info[i]['Real_str'], 
                            dirpath + '/%s_spin_correlations.inc' % overall_sector_info[i]['Real_str'] )


    #===========================================================================
    # write file for testing limits, 'testRR.f'
    #===========================================================================

    def write_testRR_3p_template_file(self, writer, dirpath, dirmadnklo, defining_process, 
                                    i, isec, jsec, ksec, lsec, K1_ct, K2_ct,K12_ct):
        replace_dict = {}
        replace_dict['isec'] = isec
        replace_dict['jsec'] = jsec
        replace_dict['ksec'] = ksec
        replace_dict['lsec'] = lsec
        mass_list = []
        mass_list = 'ZERO' #FIRST ATTEMPT TO SKIP THE IF BELOW

        limit_str = ''
    

        # Test for K1_3p (ijk)  
        # Identify cts for L1_ijk (list starts from 0)
        # L1_ijk  : 6  -> [Si, Sj, Sk, HCij, HCik, HCjk]    
        # mapping: ((isec,jsec,iref),(jsec,ksec,iref))
        if K1_ct[i][0] != 0 : #Si
            limit_str += """
c
c     soft limit for isec particle going soft
      e = [0d0,0d0,0d0,1d0,1d0] ! Si limit
      l = [0d0,0d0,0d0,0d0,0d0]
      call do_limit_RR_%d_%d_%d_%d(iunit,'Si      ',x0,e,l)
"""%(isec,jsec,ksec,lsec)
        if K1_ct[i][1] != 0 : #Sj
            limit_str += """
c
      e=[0d0,0d0,0d0,1d0,1d0] ! Sj limit
      l=[0d0,0d0,0d0,1d0,0d0]
      call do_limit_RR_%d_%d_%d_%d(iunit,'Sj      ',x0,e,l)
"""%(isec,jsec,ksec,lsec)
        if K1_ct[i][2] != 0 : #Sk to CHECK
            limit_str += """
c
      e=[1d0,1d0,0d0,0d0,0d0] ! Sk limit
      l=[1d0,0d0,0d0,0d0,0d0]
      call do_limit_RR_%d_%d_%d_%d(iunit,'Sk      ',x0,e,l)
"""%(isec,jsec,ksec,lsec)
        # Loop over sectors with final state particles only
        if isec > 2 and jsec > 2:
            if K1_ct[i][3] != 0 : # Cij
                limit_str += """
c
c     collinear limit Cij
        e=[0d0,0d0,0d0,0d0,1d0]
        l=[0d0,0d0,0d0,0d0,0d0]
      call do_limit_RR_%d_%d_%d_%d(iunit,'Cij     ',x0,e,l)
"""%(isec,jsec,ksec,lsec)
            if K1_ct[i][4] != 0 : # Cik TO CHECK:
                limit_str += """
c
c     collinear limit Cik ! todo
      e=[0d0,0d0,0d0,0d0,0d0]
      l=[0d0,0d0,0d0,0d0,0d0]  
      call do_limit_RR_%d_%d_%d_%d(iunit,'Cik     ',x0,e,l)
"""%(isec,jsec,ksec,lsec)   
            if K1_ct[i][-1] != 0 :
                limit_str += """
c
c     collinear limit Cjk ! todo
      e=[0d0,0d0,0d0,0d0,0d0]
      l=[0d0,0d0,0d0,0d0,0d0]  
      call do_limit_RR_%d_%d_%d_%d(iunit,'Cjk     ',x0,e,l)
"""%(isec,jsec,ksec,lsec)
        elif isec > 2 and jsec <= 2:
            limit_str += """Collinear limits still to be specified in sectorsRR.py """
            raise MadEvent7Error('Collinear limits still to be specified in sectorsRR.py. ')


        # Test for K2_3p(ijk) 
        # L2_ijk  : 10 -> [Sij, Sik, Sjk,
                    #                  SHCijk, SHCjik, SHCkij, 
                    #                  HCijk, 
                    #                  CijkSHCijk, CijkSHCjik, CijkSHCkij]  
        # mapping ((isec,ksec,iref),(jsec,ksec,iref))
                    
        if K2_ct[i][0] != 0: # Sij limit
            limit_str += """
c
c     double soft limit for (isec,jsec) particles going soft
c   mapping ((i,j,r),(j,k,r))
      e = [1d0,1d0,0d0,0d0,1d0] ! Sij limit
      l = [0d0,0d0,0d0,0d0,0d0]
      call do_limit_RR_%d_%d_%d_%d(iunit,'Sij     ',x0,e,l)
"""%(isec,jsec,ksec,lsec)
        if K2_ct[i][1] != 0: # Sik limit to check
            limit_str += """
    c
c     double soft limit for (isec,ksec) particles going soft
      e = [1d0,1d0,0d0,1d0,1d0] ! Sik limit
      l = [1d0,0d0,0d0,0d0,0d0]
      call do_limit_RR_%d_%d_%d_%d(iunit,'Sik     ',x0,e,l)
"""%(isec,jsec,ksec,lsec)
        if K2_ct[i][2] != 0: # Sjk limit to check
            limit_str += """
    c
c     double soft limit for (jsec,ksec) particles going soft
      e = [1d0,1d0,0d0,1d0,1d0] ! Sjk limit
      l = [1d0,0d0,0d0,1d0,0d0]
      call do_limit_RR_%d_%d_%d_%d(iunit,'Sjk     ',x0,e,l)
"""%(isec,jsec,ksec,lsec)
        if isec > 2 and jsec > 2:
            if K2_ct[i][3] != 0 : # SHCijk
                limit_str += """
c
c     single soft double collinear limit SCijk
        e=[0d0,1d0,0d0,1d0,1d0]
        l=[0d0,0d0,0d0,0d0,0d0]
      call do_limit_RR_%d_%d_%d_%d(iunit,'SCijk ',x0,e,l)
"""%(isec,jsec,ksec,lsec)
            if K2_ct[i][4] != 0 : # SHCjik
                limit_str += """
c
c     single soft double collinear limit SCjik
        e=[0d0,1d0,0d0,1d0,1d0]
        l=[0d0,0d0,0d0,1d0,0d0]
      call do_limit_RR_%d_%d_%d_%d(iunit,'SCjik ',x0,e,l)
"""%(isec,jsec,ksec,lsec)
                if K2_ct[i][5] != 0 : # SHCkij
                    limit_str += """
c
c     single soft double collinear limit SCkij
        e=[1d0,1d0,0d0,0d0,1d0]
        l=[1d0,0d0,0d0,0d0,0d0]
      call do_limit_RR_%d_%d_%d_%d(iunit,'SCkij ',x0,e,l)
"""%(isec,jsec,ksec,lsec)
                if K2_ct[i][6] != 0 : # HCijk
                    limit_str += """
c
c       double collinear limit Cijk
        e=[0d0,1d0,0d0,0d0,1d0]
        l=[0d0,0d0,0d0,0d0,0d0]
      call do_limit_RR_%d_%d_%d_%d(iunit,'Cijk ',x0,e,l)
"""%(isec,jsec,ksec,lsec)
        replace_dict['limit_str'] = limit_str
        replace_dict['NNLO_proc_str'] = str(defining_process.shell_string(schannel=True, 
                                        forbid=True, main=False, pdg_order=False, print_id = False) + '_')
#         replace_dict['mapping_str'] = mapping_str

        # write testR             
        filename = pjoin(dirpath, 'testRR_%d_%d_%d.f' %(isec,jsec,ksec) )
        file = open(pjoin(dirmadnklo,"tmp_fortran/tmp_files/testRR_template.f")).read()
        file = file % replace_dict
        writer(filename).writelines(file)

        return True
    

    def write_testRR_4p_template_file(self, writer, dirpath, dirmadnklo, defining_process, 
                                    i, isec, jsec, ksec, lsec, K1_ct, K2_ct,K12_ct):
        replace_dict = {}
        replace_dict['isec'] = isec
        replace_dict['jsec'] = jsec
        replace_dict['ksec'] = ksec
        replace_dict['lsec'] = lsec
        mass_list = []
        mass_list = 'ZERO' #FIRST ATTEMPT TO SKIP THE IF BELOW

        limit_str = ''
        
        # Test for K1_4p (ijkl)  
        # Identify cts for L1_ijkl (list starts from 0)
        # L1_ijkl  : 6 -> [Si, Sj, Sk, Sl, HCij, HCkl]  
        # mapping ((isec,ksec,lsec),(jsec,lsec,iref))
        if K1_ct[i][0] != 0 : #Si
            limit_str += """
c
c     soft limit for isec particle going soft
      e = [1d0,1d0,0d0,0d0,0d0] ! Si limit
      l = [0d0,0d0,0d0,0d0,0d0]
      call do_limit_RR_%d_%d_%d_%d(iunit,'Si      ',x0,e,l)
"""%(isec,jsec,ksec,lsec)
        if K1_ct[i][1] != 0 : #Sj
            limit_str += """
c
      e=[0d0,0d0,0d0,1d0,1d0] ! Sj limit
      l=[0d0,0d0,0d0,0d0,0d0]
      call do_limit_RR_%d_%d_%d_%d(iunit,'Sj      ',x0,e,l)
"""%(isec,jsec,ksec,lsec)
        if K1_ct[i][2] != 0 : #Sk to do
            limit_str += """
c
c     soft limit
      e=[1d0,1d0,0d0,0d0,0d0] ! Sk limit
      l=[1d0,0d0,0d0,0d0,0d0]
      call do_limit_RR_%d_%d_%d_%d(iunit,'Sk      ',x0,e,l)
"""%(isec,jsec,ksec,lsec)
        if K1_ct[i][3] != 0 : #Sl to do
            limit_str += """
c
c     soft limit
      e=[0d0,0d0,0d0,1d0,1d0] ! Sl limit
      l=[0d0,0d0,0d0,1d0,0d0]
      call do_limit_RR_%d_%d_%d_%d(iunit,'Sl      ',x0,e,l)
"""%(isec,jsec,ksec,lsec)    
        # Loop over sectors with final state particles only
        if isec > 2 and jsec > 2:
            if K1_ct[i][4] != 0 : # Cij
                limit_str += """
c
c     use mapping ((isec,jsec,lsec),(ksec,lsec,iref))  (j <--> k)
      iS1tmp = iS1
      iU2tmp = iU2
      iS1 = iU2
      iU2 = iS1tmp
c     collinear limit Cij
        e=[0d0,1d0,0d0,0d0,0d0]
        l=[0d0,0d0,0d0,0d0,0d0]
      call do_limit_RR_%d_%d_%d_%d(iunit,'Cij     ',x0,e,l)
      iS1 = iS1tmp
      iU2 = iU2tmp
"""%(isec,jsec,ksec,lsec)
                if  K1_ct[i][0] != 0:   # SiCij 
                    limit_str += """
c
c     use mapping ((isec,jsec,lsec),(ksec,lsec,iref))  (j <--> k)
      iS1tmp = iS1
      iU2tmp = iU2
      iS1 = iU2
      iU2 = iS1tmp
c     soft-collinear limit
      e=[1d0,2d0,0d0,0d0,0d0]  
      l=[0d0,0d0,0d0,0d0,0d0]
      call do_limit_RR_%d_%d_%d_%d(iunit,'SiCij      ',x0,e,l)
      iS1 = iS1tmp
      iU2 = iU2tmp
"""%(isec,jsec,ksec,lsec)
                if K1_ct[i][1] != 0 : #SjCij
                    limit_str += """
c
c     use mapping ((isec,jsec,lsec),(ksec,lsec,iref))  (j <--> k)
      iS1tmp = iS1
      iU2tmp = iU2
      iS1 = iU2
      iU2 = iS1tmp
c     soft-collinear limit
      e=[1d0,2d0,0d0,0d0,0d0]
      l=[1d0,0d0,0d0,0d0,0d0]
      call do_limit_RR_%d_%d_%d_%d(iunit,'SjCij      ',x0,e,l)
      iS1 = iS1tmp
      iU2 = iU2tmp
      """%(isec,jsec,ksec,lsec)
            if K1_ct[i][5] != 0 : # Ckl:
                limit_str += """
      
c
c     collinear limit Ckl
      e=[0d0,0d0,0d0,0d0,1d0]
      l=[0d0,0d0,0d0,0d0,1d0]  
      call do_limit_RR_%d_%d_%d_%d(iunit,'Ckl     ',x0,e,l)
"""%(isec,jsec,ksec,lsec)
                if K1_ct[i][2] != 0 : # SkCkl
                    limit_str += """
c
c     soft-collinear limit
      e=[1d0,0d0,0d0,0d0,2d0]
      l=[1d0,0d0,0d0,0d0,1d0] 
      call do_limit_RR_%d_%d_%d_%d(iunit,'SkCkl      ',x0,e,l)
"""%(isec,jsec,ksec,lsec)
                if K1_ct[i][3] != 0 : # SlCkl
                    limit_str += """
c
c     soft-collinear limit
      e=[0d0,0d0,0d0,1d0,2d0]
      l=[0d0,0d0,0d0,1d0,1d0] 
      call do_limit_RR_%d_%d_%d_%d(iunit,'SlCkl      ',x0,e,l)
"""%(isec,jsec,ksec,lsec) 
                
        elif isec > 2 and jsec <= 2:
            limit_str += """Collinear limits still to be specified in sectorsRR.py """
            raise MadEvent7Error('Collinear limits still to be specified in sectorsRR.py. ')



        replace_dict['limit_str'] = limit_str
        replace_dict['NNLO_proc_str'] = str(defining_process.shell_string(schannel=True, 
                                        forbid=True, main=False, pdg_order=False, print_id = False) + '_')
#         replace_dict['mapping_str'] = mapping_str

        # write testR             
        filename = pjoin(dirpath, 'testRR_%d_%d_%d_%d.f' %(isec,jsec,ksec,lsec) )
        file = open(pjoin(dirmadnklo,"tmp_fortran/tmp_files/testRR_template.f")).read()
        file = file % replace_dict
        writer(filename).writelines(file)

        return True
    

    #===========================================================================
    # write driver_isec_jsec for real subprocess directory
    #===========================================================================

    def write_driver_npt_template(self, writer, dirpath, dirmadnklo, i , isec, jsec, ksec, lsec, UBgraphs):
        
        replace_dict = {}
        replace_dict['isec'] = isec
        replace_dict['jsec'] = jsec
        replace_dict['ksec'] = ksec
        replace_dict['lsec'] = lsec
        #replace_dict['UBgraphs'] = UBgraphs

        # write driver
        if(lsec != 0):
            filename = pjoin(dirpath, 'driver_RR_%d_%d_%d_%d.f' % (isec, jsec,ksec,lsec))
        else:
            filename = pjoin(dirpath, 'driver_RR_%d_%d_%d.f' % (isec, jsec,ksec))
        file = open(pjoin(dirmadnklo,"tmp_fortran/tmp_files/driver_npt_template.f")).read()
        file = file % replace_dict
        writer(filename).writelines(file)

        return True

    
    def write_makefile_RR_file(self, writer, dirpath, dirmadnklo, defining_process, overall_sector_info):

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
    # write driver_isec_jsec for real subprocess directory
    #===========================================================================

    def write_driver_npt_template(self, writer, dirpath, dirmadnklo, i , isec, jsec, ksec, lsec, UBgraphs):
        
        replace_dict = {}
        replace_dict['isec'] = isec
        replace_dict['jsec'] = jsec
        replace_dict['ksec'] = ksec
        replace_dict['lsec'] = lsec
        #replace_dict['UBgraphs'] = UBgraphs

        # write driver
        if(lsec != 0):
            filename = pjoin(dirpath, 'driver_RR_%d_%d_%d_%d.f' % (isec, jsec,ksec,lsec))
        else:
            filename = pjoin(dirpath, 'driver_RR_%d_%d_%d.f' % (isec, jsec,ksec))
        file = open(pjoin(dirmadnklo,"tmp_fortran/tmp_files/driver_npt_template.f")).read()
        file = file % replace_dict
        writer(filename).writelines(file)

        return True

    
    def write_makefile_RR_file(self, writer, dirpath, dirmadnklo, defining_process, overall_sector_info):

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
