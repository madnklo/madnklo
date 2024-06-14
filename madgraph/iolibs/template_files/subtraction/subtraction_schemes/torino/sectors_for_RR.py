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
        necessary_3p_ct1_list = [0] * (6*len(all_3p_sectors))
        necessary_3p_ct1 = [0] * (6*len(all_3p_sectors))
        necessary_3p_ct2_list = [0] * (10*len(all_3p_sectors))
        necessary_3p_ct2 = [0] * (10*len(all_3p_sectors))
        necessary_3p_ct12_list = [0] * (24*len(all_3p_sectors))
        necessary_3p_ct12 = [0] * (24*len(all_3p_sectors))
        i = 0
        for s in all_3p_sectors:
            s['sector'].all_3p_sector_list = all_3p_sector_list
            s['sector'].all_3p_sector_id_list = all_3p_sector_id_list

            if counterterms is not None:
                s['counterterms'] = []
                for i_ct, ct in enumerate(counterterms):
                    print('ict + ct : ' + str(i_ct) + ' ' + str(ct))
                    current = ct.nodes[0].current
                    n_subs = current.get('singular_structure').substructures
                    singular_structure = current.get('singular_structure').substructures[0]
                    all_legs = singular_structure.get_all_legs()
                    leg_numbers = s['sector'].leg_numbers
                    ileg = leg_numbers[0]
                    jleg = leg_numbers[1]
                    kleg = leg_numbers[2]
                    print(str(ileg) + ' ' + str(jleg) + ' ' + str(kleg))
                    #print(str(all_3p_sector_id_list))

                    # safety check
                    if not sorted([l.n for l in current.get('singular_structure').get_all_legs()]) == sorted([ileg,jleg,kleg]) :
                        continue

                    #print('singular_structure.name() : ' + str(singular_structure.name()))
                    #if singular_structure.substructures:
                    #    print('singular_structure        : ' + str(singular_structure.substructures[0]))

                    if len(current.get('singular_structure').substructures) > 1 :
                        print('# of subs : ' + str(len(current.get('singular_structure').substructures)))
                        for i in range(0,len(current.get('singular_structure').substructures)):
                            #print(i)
                            subs = current.get('singular_structure').substructures[i-1]
                            print('subs ' + ' : ' + str(subs))

                            n = len(subs.substructures) 
                            for j in range(1,n) :
                                print('subs.sub : ' + str(subs.substructures[j-1]))
                    else :
                        print('singular_structure  : ' + str(singular_structure))
                        n = len(singular_structure.substructures) 
                        if n > 0 :
                            for k in range(0,n) :
                                print('subs.sub1 : ' + str(singular_structure.substructures[k-1]))
                                m = len(singular_structure.substructures[k-1].substructures) 
                                if m > 0 :
                                    print('subs.sub2 : ' + str(singular_structure.substructures[k-1].substructures[0]))
                                    

                    #    print('subtructures              : ' + str(current.get('singular_structure').substructures[0]))
                    #    print('subtructures              : ' + str(current.get('singular_structure').substructures[1]))
                    #print('all legs                  : ' + str(current.get('singular_structure').get_all_legs()))

                    # Identify cts for L1_ijk
                    # L1_ijk  : 6  -> [Si, Sj, Sk, HCij, HCik, HCjk]

                    # Si
                    if s['sector'].id[0] == 21 :
                        s['counterterms'].append(i_ct)
                        necessary_3p_ct1_list[1] = 1
                        necessary_3p_ct1[1] = ct
                    # Sj
                    if s['sector'].id[1] == 21 :
                        s['counterterms'].append(i_ct)
                        necessary_3p_ct1_list[2] = 1
                        necessary_3p_ct1[2] = ct
                    # Sk
                    if s['sector'].id[2] == 21 :
                        s['counterterms'].append(i_ct)
                        necessary_3p_ct1_list[3] = 1
                        necessary_3p_ct1[3] = ct
                    
                    if singular_structure.name()=='C' and len(all_legs)==2:
                        if not singular_structure.substructures:
                            # Cij
                            if sorted([l.n for l in all_legs]) == (sorted([ileg,jleg])):
                                s['counterterms'].append(i_ct)
                                necessary_3p_ct1_list[4] = 1
                                necessary_3p_ct1[4] = ct
                            # Cik
                            if sorted([l.n for l in all_legs]) == (sorted([ileg,kleg])):
                                s['counterterms'].append(i_ct)
                                necessary_3p_ct1_list[5] = 1
                                necessary_3p_ct1[5] = ct
                            # Cjk
                            if sorted([l.n for l in all_legs]) == (sorted([jleg,kleg])):
                                s['counterterms'].append(i_ct)
                                necessary_3p_ct1_list[6] = 1
                                necessary_3p_ct1[6] = ct


                    # Identify cts for L2_ijk
                    # L2_ijk  : 10 -> [Sij, Sik, Sjk,
                    #                  SHCijk, SHCjik, SHCkij, 
                    #                  HCijk, 
                    #                  CijkSHCijk, CijkSHCjik, CijkSHCkij]  
                                
                    if n_subs == 1 and len(singular_structure.substructures) == 0 :

                        if singular_structure.name()=='S' and len(all_legs)==2:
                            if not singular_structure.substructures:
                                # Sij
                                if sorted([l.n for l in all_legs]) == (sorted([ileg,jleg])):
                                    s['counterterms'].append(i_ct)
                                    necessary_3p_ct2_list[1] = 1
                                    necessary_3p_ct2[1] = ct
                                # Sik
                                if sorted([l.n for l in all_legs]) == (sorted([ileg,kleg])):
                                    s['counterterms'].append(i_ct)
                                    necessary_3p_ct2_list[2] = 1
                                    necessary_3p_ct2[2] = ct
                                # Sjk
                                if sorted([l.n for l in all_legs]) == (sorted([jleg,kleg])):
                                    s['counterterms'].append(i_ct)
                                    necessary_3p_ct2_list[3] = 1
                                    necessary_3p_ct2[3] = ct

                        if singular_structure.name()=='C' and len(all_legs)==3:
                            # Cijk
                            if sorted([l.n for l in all_legs]) == (sorted([ileg,jleg,kleg])):
                                s['counterterms'].append(i_ct)
                                necessary_3p_ct2_list[7] = 1
                                necessary_3p_ct2[7] = ct

                    if n_subs == 2 : 
                        # here singular_structure = coll_sub -> C(i,j)

                        # SCijk 
                        if sorted([l.n for l in all_legs]) == (sorted([jleg,kleg])) :
                            s['counterterms'].append(i_ct)
                            necessary_3p_ct2_list[4] = 1
                            necessary_3p_ct2[4] = ct
                        # SCjik
                        if sorted([l.n for l in all_legs]) == (sorted([ileg,kleg])) :
                            s['counterterms'].append(i_ct)
                            necessary_3p_ct2_list[5] = 1
                            necessary_3p_ct2[5] = ct
                        # SCkij
                        if sorted([l.n for l in all_legs]) == (sorted([ileg,jleg])) :
                            s['counterterms'].append(i_ct)
                            necessary_3p_ct2_list[6] = 1
                            necessary_3p_ct2[6] = ct

                    if n_subs == 1 and len(singular_structure.substructures) == 2 :
                        coll_subsub = singular_structure.substructures[0]

                        # CijkSCijk
                        if sorted([l.n for l in coll_subsub.get_all_legs()]) == (sorted([jleg,kleg])) :
                            s['counterterms'].append(i_ct)
                            necessary_3p_ct2_list[8] = 1
                            necessary_3p_ct2[8] = ct
                        # CijkSCjik
                        if sorted([l.n for l in coll_subsub.get_all_legs()]) == (sorted([ileg,kleg])) :
                            s['counterterms'].append(i_ct)
                            necessary_3p_ct2_list[9] = 1
                            necessary_3p_ct2[9] = ct
                        # CijkSCkij
                        if sorted([l.n for l in coll_subsub.get_all_legs()]) == (sorted([ileg,jleg])) :
                            s['counterterms'].append(i_ct)
                            necessary_3p_ct2_list[10] = 1
                            necessary_3p_ct2[10] = ct

                    # NOTE: missing SijCijk, SijSCkij 

                    # Identify cts for L12_ijk
                    # L12_ijk : 24 -> [Si Sij, Si Sik, Si SHCijk, Si HCijk',
                    #                  Sj Sij, Sj Sjk, Sj SHCjik, Sj HCijk', 
                    #                  Sk Sik, Sk Sjk, Sk SHCkij, Sk HCijk',     12              
                    #                  HCij Sij, HCij SCkij, HCij Cijk(1-Sij), HCij CijkSCkij, 
                    #                  HCik Sik, HCik SCjik, HCik Cijk(1-Sik), HCik CijkSCjik, 
                    #                  HCjk Sjk, HCjk SCijk, HCjk Cijk(1-Sjk), HCjk CijkSCijk]  12                    

                    if s['sector'].id[0] == 21 and len(all_legs)>=2:
                        # Si Sij
                        necessary_3p_ct12_list[1] = necessary_3p_ct2_list[1]*1
                        # Si Sik
                        necessary_3p_ct12_list[2] = necessary_3p_ct2_list[2]*1
                        # Si SHCijk = Si SCijk(1-Sij-Sik)
                        necessary_3p_ct12_list[3] = necessary_3p_ct2_list[4]*1
                        # Si HCijk' = Si Cijk(1-Sij-Sik)
                        necessary_3p_ct12_list[4] = necessary_3p_ct2_list[7]*1

                    if s['sector'].id[1] == 21 and len(all_legs)>=2:
                        # Sj Sij
                        necessary_3p_ct12_list[5] = necessary_3p_ct2_list[1]*1
                        # Sj Sjk
                        necessary_3p_ct12_list[6] = necessary_3p_ct2_list[3]*1
                        # Sj SHCjik = Sj SCjik(1-Sij-Sjk)
                        necessary_3p_ct12_list[7] = necessary_3p_ct2_list[5]*1
                        # Sj HCijk' = Sj Cijk(1-Sij-Sjk)
                        necessary_3p_ct12_list[8] = necessary_3p_ct2_list[7]*1

                    if s['sector'].id[2] == 21 and len(all_legs)>=2:
                        # Sk Sik
                        necessary_3p_ct12_list[9] = necessary_3p_ct2_list[2]*1
                        # Sk Sjk
                        necessary_3p_ct12_list[10] = necessary_3p_ct2_list[3]*1
                        # Sk SHCkij = Sk SCkij(1-Sik-Sjk)
                        necessary_3p_ct12_list[11] = necessary_3p_ct2_list[6]*1
                        # Sk HCijk' = Sk Cijk(1-Sik-Sjk)
                        necessary_3p_ct12_list[12] = necessary_3p_ct2_list[7]*1

                    if n_subs == 1 and len(singular_structure.substructures) == 1 :
                        if len(all_legs)==2:
                            # HCij Sij
                            if sorted([l.n for l in all_legs]) == (sorted([ileg,jleg])):
                                necessary_3p_ct12_list[13] = necessary_3p_ct2_list[1]*1
                            # HCik Sik
                            if sorted([l.n for l in all_legs]) == (sorted([ileg,jleg])):
                                necessary_3p_ct12_list[17] = necessary_3p_ct2_list[2]*1
                            # HCjk Sjk
                            if sorted([l.n for l in all_legs]) == (sorted([ileg,jleg])):
                                necessary_3p_ct12_list[21] = necessary_3p_ct2_list[3]*1

                        # HCij SCkij
                        necessary_3p_ct12_list[14] = necessary_3p_ct2_list[6]*1
                        # HCik SCjik
                        necessary_3p_ct12_list[18] = necessary_3p_ct2_list[5]*1
                        # HCjk SCijk
                        necessary_3p_ct12_list[22] = necessary_3p_ct2_list[4]*1
                        # HCij Cijk
                        necessary_3p_ct12_list[15] = necessary_3p_ct2_list[7]*1
                        # HCik Cijk
                        necessary_3p_ct12_list[19] = necessary_3p_ct2_list[7]*1
                        # HCjk Cijk 
                        necessary_3p_ct12_list[23] = necessary_3p_ct2_list[7]*1
                        # HCij CijkSCkij
                        necessary_3p_ct12_list[16] = necessary_3p_ct2_list[10]*1
                        # HCik CijkSCjik
                        necessary_3p_ct12_list[20] = necessary_3p_ct2_list[9]*1
                        # HCjk CijkSCijk
                        necessary_3p_ct12_list[24] = necessary_3p_ct2_list[8]*1

        ######## 4p #########
        all_4p_local_counterterms_list = []
        necessary_4p_ct1_list = [0] * (6*len(all_4p_sectors))
        necessary_4p_ct1 = [0] * (6*len(all_4p_sectors))
        necessary_4p_ct2_list = [0] * (10*len(all_4p_sectors))
        necessary_4p_ct2 = [0] * (13*len(all_4p_sectors))
        necessary_4p_ct12_list = [0] * (24*len(all_4p_sectors))
        necessary_4p_ct12 = [0] * (18*len(all_4p_sectors))
        for s in all_4p_sectors:
            s['sector'].all_4p_sector_list = all_4p_sector_list
            s['sector'].all_4p_sector_id_list = all_4p_sector_id_list

            if counterterms is not None:
                s['counterterms'] = []
                for i_ct, ct in enumerate(counterterms):
                    print('ict + ct : ' + str(i_ct) + ' ' + str(ct))
                    current = ct.nodes[0].current
                    n_subs = current.get('singular_structure').substructures
                    singular_structure = current.get('singular_structure').substructures[0]
                    all_legs = singular_structure.get_all_legs()
                    leg_numbers = s['sector'].leg_numbers
                    ileg = leg_numbers[0]
                    jleg = leg_numbers[1]
                    kleg = leg_numbers[2]
                    lleg = leg_numbers[3]
                    print(str(ileg) + ' ' + str(jleg) + ' ' + str(kleg) + ' ' + str(lleg))


                    # L1_ijkl  : 6 -> [Si, Sj, Sk, Sl, HCij, HCkl]

                    # Si
                    if s['sector'].id[0] == 21 :
                        s['counterterms'].append(i_ct)
                        necessary_4p_ct1_list[1] = 1
                        necessary_4p_ct1[1] = ct
                    # Sj
                    if s['sector'].id[1] == 21 :
                        s['counterterms'].append(i_ct)
                        necessary_4p_ct1_list[2] = 1
                        necessary_4p_ct1[2] = ct
                    # Sk
                    if s['sector'].id[2] == 21 :
                        s['counterterms'].append(i_ct)
                        necessary_4p_ct1_list[3] = 1
                        necessary_4p_ct1[3] = ct
                    # Sl
                    if s['sector'].id[3] == 21 :
                        s['counterterms'].append(i_ct)
                        necessary_4p_ct1_list[4] = 1
                        necessary_4p_ct1[4] = ct
                    
                    if singular_structure.name()=='C' and len(all_legs)==2:
                        if not singular_structure.substructures:
                            # Cij
                            if sorted([l.n for l in all_legs]) == (sorted([ileg,jleg])):
                                s['counterterms'].append(i_ct)
                                necessary_4p_ct1_list[5] = 1
                                necessary_4p_ct1[5] = ct
                            # Ckl
                            if sorted([l.n for l in all_legs]) == (sorted([kleg,lleg])):
                                s['counterterms'].append(i_ct)
                                necessary_4p_ct1_list[6] = 1
                                necessary_4p_ct1[6] = ct

                    # L2_ijkl  : 13  ->  [Sik, Sil, Sjk, Sjl, 
                    #                     SHCikl, SHCjkl, SHCkij, SHClij, 
                    #                     Cijkl(1+Sik+Sil+Sjk+Sjl), 
                    #                     CijklSCikl, CijklSCjkl, CijklSCkij, CijklSClij]  

                    if n_subs == 1 and len(singular_structure.substructures) == 0 :

                        if singular_structure.name()=='S' and len(all_legs)==2:
                            if not singular_structure.substructures:
                                # Sik
                                if sorted([l.n for l in all_legs]) == (sorted([ileg,kleg])):
                                    s['counterterms'].append(i_ct)
                                    necessary_4p_ct2_list[1] = 1
                                    necessary_4p_ct2[1] = ct
                                # Sil
                                if sorted([l.n for l in all_legs]) == (sorted([ileg,lleg])):
                                    s['counterterms'].append(i_ct)
                                    necessary_4p_ct2_list[2] = 1
                                    necessary_4p_ct2[2] = ct
                                # Sjk
                                if sorted([l.n for l in all_legs]) == (sorted([jleg,kleg])):
                                    s['counterterms'].append(i_ct)
                                    necessary_4p_ct2_list[3] = 1
                                    necessary_4p_ct2[3] = ct
                                # Sjl
                                if sorted([l.n for l in all_legs]) == (sorted([jleg,lleg])):
                                    s['counterterms'].append(i_ct)
                                    necessary_4p_ct2_list[4] = 1
                                    necessary_4p_ct2[4] = ct

                        if singular_structure.name()=='C' and len(all_legs)==4:
                            # Cijkl
                            if sorted([l.n for l in all_legs]) == (sorted([ileg,jleg,kleg,lleg])):
                                s['counterterms'].append(i_ct)
                                necessary_4p_ct2_list[7] = 1
                                necessary_4p_ct2[7] = ct

                    if n_subs == 2 : 
                        # here singular_structure = coll_sub -> C(i,j)

                        # CijklSCikl
                        if sorted([l.n for l in all_legs]) == (sorted([jleg,kleg])) :
                            s['counterterms'].append(i_ct)
                            necessary_3p_ct2_list[4] = 1
                            necessary_3p_ct2[4] = ct
                        # SCjik
                        if sorted([l.n for l in all_legs]) == (sorted([ileg,kleg])) :
                            s['counterterms'].append(i_ct)
                            necessary_3p_ct2_list[5] = 1
                            necessary_3p_ct2[5] = ct
                        # SCkij
                        if sorted([l.n for l in all_legs]) == (sorted([ileg,jleg])) :
                            s['counterterms'].append(i_ct)
                            necessary_3p_ct2_list[6] = 1
                            necessary_3p_ct2[6] = ct

                    if n_subs == 1 and len(singular_structure.substructures) == 2 :
                        coll_subsub = singular_structure.substructures[0]

                        # CijkSCijk
                        if sorted([l.n for l in coll_subsub.get_all_legs()]) == (sorted([jleg,kleg])) :
                            s['counterterms'].append(i_ct)
                            necessary_3p_ct2_list[8] = 1
                            necessary_3p_ct2[8] = ct
                        # CijkSCjik
                        if sorted([l.n for l in coll_subsub.get_all_legs()]) == (sorted([ileg,kleg])) :
                            s['counterterms'].append(i_ct)
                            necessary_3p_ct2_list[9] = 1
                            necessary_3p_ct2[9] = ct
                        # CijkSCkij
                        if sorted([l.n for l in coll_subsub.get_all_legs()]) == (sorted([ileg,jleg])) :
                            s['counterterms'].append(i_ct)
                            necessary_3p_ct2_list[10] = 1
                            necessary_3p_ct2[10] = ct

                    # L12_ijkl : 18  ->  [Si Sik, Sk Sik, Si Sil, Sl Sil, Sj Sjk, Sk Sjk, Sj Sjl, Sl Sjl,
                    #                     Si SHCikl, Sj SHCjkl, Sk SHCkij, Sl SHClij,
                    #                     HCij SCkij, HCij SClij, HCkl SCikl, HCkl SCjkl,
                    #                     HCij Cijkl, HCkl Cijkl,
                    #                     HCij CijklSCkij, HCij CijklSClij, HCkl CijklSCikl, HCkl CijklSCjkl] 


                    # TODO: finish and check 4p sectors
                            
        # Find the corresponding integrated counterterms (None for RR)



######### Set writer
        writer = writers.FortranWriter

        # Point to the right process directory
#        dirmadnklo=os.getcwd()
#        dirpath = pjoin(dirmadnklo,glob.glob("%s/NNLO_RR_x_RR_*" % interface.user_dir_name[0])[0])
#        dirpath = pjoin(dirpath, 'SubProcesses', \
#                       "P%s" % defining_process.shell_string())
#        sys.path.append(pjoin(dirmadnklo,"%s/SubProcesses" % interface.user_dir_name[0]))

######### Write all_sector_list.inc
        #self.write_all_sector_list_include(writers.FortranWriter, dirpath, all_sector_list)
            
######### Write NNLO_K_isec_jsec_ksec.f (3-particle sector)
        
        # Set replace_dict for NLO_K_isec_jsec.f
        # replace_dict_ct = {}
        # replace_dict_limits = {}
        # for i in range(0,len(all_3p_sector_list)):
        #     list_M2 = []
        #     list_str_defHC = []
        #     mapping = []
        #     mapping_str = ''
        #     isec = all_3p_sector_list[i][0]
        #     jsec = all_3p_sector_list[i][1]
        #     ksec = all_3p_sector_list[i][2]
        #     id_isec = all_3p_sector_id_list[i][0]
        #     id_jsec = all_3p_sector_id_list[i][1]
        #     id_ksec = all_3p_sector_id_list[i][2]
        #     # Extract the reference particle leg from recoiler_function.py
        #     iref = all_sector_recoilers[i]
        #     replace_dict_ct['iref'] = iref
        #     # Check isec != jsec
        #     if isec == jsec:
        #         raise MadEvent7Error('Wrong sector indices %d,%d!' % (isec,jsec))
        #     replace_dict_ct['isec'] = isec
        #     replace_dict_ct['jsec'] = jsec
        #     replace_dict_ct['ksec'] = ksec
        #     replace_dict_limits['isec'] = isec
        #     replace_dict_limits['jsec'] = jsec
        #     replace_dict_limits['ksec'] = ksec
        #     replace_dict_limits['proc_prefix_real'] = str(defining_process.shell_string(schannel=True, 
        #                                 forbid=True, main=False, pdg_order=False, print_id = False))
            
        #     # 5 strutture:  ggg  ggq  gqqx  qqxq'  qqxq     
        #     # K1 (4)
        #     # Si
        #     replace_dict_limits['proc_prefix_S_g'] = 'dummy'
        #     # HCij
        #     replace_dict_limits['proc_prefix_HC_gg'] = 'dummy'
        #     replace_dict_limits['proc_prefix_HC_gq'] = 'dummy'
        #     replace_dict_limits['proc_prefix_HC_qqx'] = 'dummy'

        #     # K2 (15)
        #     # Sij
        #     replace_dict_limits['proc_prefix_SS_gg'] = 'dummy'
        #     replace_dict_limits['proc_prefix_SS_qqx'] = 'dummy'
        #     # SCijk
        #     replace_dict_limits['proc_prefix_SC_ggg'] = 'dummy'
        #     replace_dict_limits['proc_prefix_SC_ggq'] = 'dummy'
        #     replace_dict_limits['proc_prefix_SC_gqqx'] = 'dummy'
        #     # Sij SCijk
        #     replace_dict_limits['proc_prefix_SC_SS_ggg'] = 'dummy'
        #     replace_dict_limits['proc_prefix_SC_SS_ggq'] = 'dummy'
        #     # Cijk
        #     replace_dict_limits['proc_prefix_HC_qqxq'] = 'dummy'
        #     replace_dict_limits['proc_prefix_HC_qqxqp'] = 'dummy'
        #     replace_dict_limits['proc_prefix_HC_gqqx'] = 'dummy'
        #     replace_dict_limits['proc_prefix_HC_ggq'] = 'dummy'
        #     replace_dict_limits['proc_prefix_HC_ggg'] = 'dummy'
        #     # CijkSHCijk
        #     replace_dict_limits['proc_prefix_CSHC_ggg'] = 'dummy'
        #     replace_dict_limits['proc_prefix_CSHC_ggq'] = 'dummy'
        #     replace_dict_limits['proc_prefix_CSHC_gqqx'] = 'dummy'

        #     # K12
        #     # SiSij
        #     replace_dict_limits['proc_prefix_S_SS_gg'] = 'dummy'
        #     # SiSHCijk
        #     replace_dict_limits['proc_prefix_S_SC_ggg'] = 'dummy'
        #     replace_dict_limits['proc_prefix_S_SC_ggq'] = 'dummy'
        #     replace_dict_limits['proc_prefix_S_SC_gqqx'] = 'dummy'
        #     # SiHC(s)ijk  where  HC(s)ijk = Cijk(1-Sij-Sik)(1-SCijk) ???
        #     replace_dict_limits['proc_prefix_S_HCs_gqqx'] = 'dummy'
        #     replace_dict_limits['proc_prefix_S_HCs_ggq'] = 'dummy'
        #     replace_dict_limits['proc_prefix_S_HCs_ggg'] = 'dummy'
        #     # HCijSij
        #     replace_dict_limits['proc_prefix_HC_SS_gg'] = 'dummy'
        #     replace_dict_limits['proc_prefix_HC_SS_qqx'] = 'dummy'
        #     # HCijSCkij
        #     replace_dict_limits['proc_prefix_HC_SC_ggg'] = 'dummy'
        #     replace_dict_limits['proc_prefix_HC_SC_ggq'] = 'dummy'
        #     replace_dict_limits['proc_prefix_HC_SC_gqqx'] = 'dummy'
        #     # HCijHC(c)ijk where HC(c)ijk = Cijk(1-Sij-SCkij)
        #     replace_dict_limits['proc_prefix_HC_HC_qqxq'] = 'dummy'
        #     replace_dict_limits['proc_prefix_HC_HC_qqxqp'] = 'dummy'
        #     replace_dict_limits['proc_prefix_HC_HC_gqqx'] = 'dummy'
        #     replace_dict_limits['proc_prefix_HC_HC_ggq'] = 'dummy'
        #     replace_dict_limits['proc_prefix_HC_HC_ggg'] = 'dummy'

    
        return #all_sectors
