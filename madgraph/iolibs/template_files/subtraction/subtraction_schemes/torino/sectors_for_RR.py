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

        print('INTO RR SECTOR')
        print(defining_process.shell_string())
        print(contrib_definition.get_shell_name())
        print(contrib_definition.process_definition.get('id'))

        all_sectors = []
        all_sector_legs = []
        all_sector_id_legs = []
        all_sector_recoilers = []
        all_3p_sector_recoilers = []

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

        # at NNLO: from i,j,k to ijk
        #      3) [g0 g1 g1 g2], (3)
        #         [g0 g1 g1 q], [g0 g1 q g1], [q <-> bq],  (2)
        #         [g q q bq], [g q bq q], (1)  (gb: [q g g bq] cannot appear, gluon always first!)
        #         [bq q1 q1 q'], [bq q1 q' q1] + [q' <-> bq'], [bq q1 q1 q2], [bq q1 q2 q1] + [q <-> bq] (0)

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
                        #a_sector['recoiler'] = None
                        a_sector['recoiler'] = recoiler_function.get_recoiler(defining_process,(i.get('number'),j.get('number')))
                        all_sector_recoilers.append(a_sector['recoiler'].get('number'))
                        all_sector_legs.append(i.get('number'))
                        all_sector_legs.append(j.get('number'))
                        # keep track of the masses
                        a_sector['sector'].masses = (model.get('particle_dict')[i.get('id')]['mass'],
                                                     model.get('particle_dict')[j.get('id')]['mass'])

                        # keep track of the particles' identity
                        a_sector['sector'].id = (i.get('id'), j.get('id'))
                        all_sector_id_legs.append(i.get('id'))
                        all_sector_id_legs.append(j.get('id'))
                        all_sectors.append(a_sector)

                    for k, col_k in zip(leglist, colorlist):
                        if k.get('number') == i.get('number') or k.get('number') == j.get('number'):
                            continue
                        if col_k == 1:
                            continue
                        if not k.get('state'):
                            continue

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
                            new_leglist = copy.copy(leglist_ij)
                            new_leglist[min([leglist_ij.index(ij), leglist_ij.index(k)])] = ijk
                            new_leglist.pop(max([leglist_ij.index(ij), leglist_ij.index(k)]))
                            new_process['legs'] = new_leglist
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
                                #a_sector['recoiler'] = None
                                a_sector['recoiler'] = recoiler_function.get_recoiler(defining_process,(i.get('number'),j.get('number')))
                                all_sector_recoilers.append(a_sector['recoiler'].get('number'))
                                all_sector_legs.append(ij.get('number'))
                                all_sector_legs.append(k.get('number'))
                                # keep track of the masses
                                a_sector['sector'].masses = (model.get('particle_dict')[ij.get('id')]['mass'],
                                                     model.get('particle_dict')[k.get('id')]['mass'])

                                # keep track of the particles' identity
                                a_sector['sector'].id = (ij.get('id'), k.get('id'))
                                all_sector_id_legs.append(ij.get('id'))
                                all_sector_id_legs.append(k.get('id'))

                                all_sectors.append(a_sector)


                        tmp_sector = [i.get('number'),j.get('number'),k.get('number')]
                        tmp_sector_id = [i.get('id'),j.get('id'),k.get('id')]

                        if len(threep_sectors) == 0 or tmp_sector not in combs:
                            ### ijjk topology ###
                            a_3p_sector = {
                            'sector': None,
                            'counterterms': None,
                            'integrated_counterterms': None,
                            'recoiler' : None
                            }
                            a_3p_sector['sector'] = sectors.Sector(leg_numbers=(i.get('number'),j.get('number'),j.get('number'),k.get('number')))
                            a_sector['recoiler'] = recoiler_function.get_recoiler(defining_process,(i.get('number'),j.get('number'),k.get('number')))
                            all_3p_sector_recoilers.append(a_sector['recoiler'].get('number'))
                            # keep track of the masses
                            #a_sector['sector'].masses = (model.get('particle_dict')[i.get('id')]['mass'],
                            #                         model.get('particle_dict')[j.get('id')]['mass'])

                            # keep track of the particles' identity
                            a_3p_sector['sector'].id = (i.get('id'),j.get('id'),j.get('id'),k.get('id'))
                            all_3p_sectors.append(a_3p_sector)
                            logger.info('NNLO sector found, legs %d, %d, %d, %d' % a_3p_sector['sector'].leg_numbers)

                            ### ijkj topology ###
                            a_3p_sector = {
                            'sector': None,
                            'counterterms': None,
                            'integrated_counterterms': None,
                            'recoiler' : None
                            }
                            a_3p_sector['sector'] = sectors.Sector(leg_numbers=(i.get('number'),j.get('number'),k.get('number'),j.get('number')))
                            a_sector['recoiler'] = recoiler_function.get_recoiler(defining_process,(i.get('number'),j.get('number'),k.get('number')))
                            all_3p_sector_recoilers.append(a_sector['recoiler'].get('number'))
                            # keep track of the masses
                            #a_sector['sector'].masses = (model.get('particle_dict')[i.get('id')]['mass'],
                            #                         model.get('particle_dict')[j.get('id')]['mass'])

                            # keep track of the particles' identity
                            a_3p_sector['sector'].id = (i.get('id'),j.get('id'),k.get('id'),j.get('id'))
                            all_3p_sectors.append(a_3p_sector)
                            logger.info('NNLO sector found, legs %d, %d, %d, %d' % a_3p_sector['sector'].leg_numbers)

                            threep_sectors.append([i.get('number'),j.get('number'),j.get('number'),k.get('number')])  #ijjk
                            threep_sectors_id.append([i.get('id'),j.get('id'),j.get('id'),k.get('id')])
                            threep_sectors.append([i.get('number'),j.get('number'),k.get('number'),j.get('number')])  #ijkj
                            threep_sectors_id.append([i.get('id'),j.get('id'),k.get('id'),j.get('id')])

                            combs.append([i.get('number'),j.get('number'),k.get('number')])
                            #combs.append([i.get('number'),j.get('number'),j.get('number'),k.get('number')])     #ijjk in combs
                            #combs.append([i.get('number'),j.get('number'),k.get('number'),j.get('number')])     #ijkj in combs

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
                        all_sector_legs.append(i.get('number'))
                        all_sector_legs.append(j.get('number'))
                        # keep track of the masses
                        a_sector['sector'].masses = (model.get('particle_dict')[i.get('id')]['mass'],
                                                     model.get('particle_dict')[j.get('id')]['mass'])

                        # keep track of the particles' identity
                        a_sector['sector'].id = (i.get('id'), j.get('id'))
                        all_sector_id_legs.append(i.get('id'))
                        all_sector_id_legs.append(j.get('id'))

                        all_sectors.append(a_sector)


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
                        if k['id'] != 21 and l['id'] == 21 and l['state']:  # [g q] [q g]
                            continue

                        # if both k and l are gluons, then keep just the case in which k (number) < l (number)
                        if k['id'] == 21 and l['id'] == 21 and l['state']: # [g1 g2] [g2 g1]
                            if l.get('number') < k.get('number') :
                                continue

                        # if j and i are quarks and antiquark in the final state, let j be the quark
                        #   this is needed in order to comply with the fct combine_ij inside fks_common
                        if k['id'] == -l['id'] and l['state']: # [q qb] [qb q]
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

                                all_sectors.append(a_sector)


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

                            # keep track of the particles' identity
                            a_4p_sector['sector'].id = (i.get('id'), j.get('id'),k.get('id'),l.get('id'))
                            all_4p_sectors.append(a_4p_sector)
                            logger.info('NNLO sector found, legs %d, %d, %d, %d' % a_4p_sector['sector'].leg_numbers)

                            fourp_sectors.append(tmp_sector)
                            fourp_sectors_id.append(tmp_sector_id)
                            combs.append([i.get('number'),j.get('number'),k.get('number'),l.get('number')])     #ijkl

                        elif tmp_sector in combs:
                            continue

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
        label = ''
        all_3p_K1_ct = []
        all_3p_K2_ct = []
        all_3p_K12_ct = []
        uB_all_3p_K1_ct = []
        uB_all_3p_K2_ct = []
        for s in all_3p_sectors:
            s['sector'].all_3p_sector_list = all_3p_sector_list
            s['sector'].all_3p_sector_id_list = all_3p_sector_id_list
            leg_numbers = s['sector'].leg_numbers

            if counterterms is not None:
                s['counterterms'] = []
                necessary_3p_ct1_list = [0] * (3)
                necessary_3p_ct1 = [0] * (3)
                necessary_3p_ct12_list = [0] * (13) # same number of CTs in K12 for ijjk&ijkj
                if (leg_numbers[1] == leg_numbers[2]):  #ijjk
                    label = 'ijjk'
                    ileg = leg_numbers[0]
                    jleg = leg_numbers[1]
                    kleg = leg_numbers[3]
                    i_id = s['sector'].id[0]
                    j_id = s['sector'].id[1]
                    k_id = s['sector'].id[3]
                    necessary_3p_ct2_list = [0] * (7)
                    necessary_3p_ct2 = [0] * (7)
                elif (leg_numbers[1] == leg_numbers[3]):  #ijkj
                    label = 'ijkj'
                    ileg = leg_numbers[0]
                    jleg = leg_numbers[1]
                    kleg = leg_numbers[2]
                    i_id = s['sector'].id[0]
                    j_id = s['sector'].id[1]
                    k_id = s['sector'].id[2]
                    necessary_3p_ct2_list = [0] * (11)
                    necessary_3p_ct2 = [0] * (11)

                print('****** NEW SECTOR ******')
                print(str(s['sector'].leg_numbers[0]) + ' ' + str(s['sector'].leg_numbers[1]) + ' ' + str(s['sector'].leg_numbers[2]) + ' ' + str(s['sector'].leg_numbers[3]))
                print(str(s['sector'].id[0]) + ' ' + str(s['sector'].id[1]) + ' ' + str(s['sector'].id[2]) + ' ' + str(s['sector'].id[3]))

                for i_ct, ct in enumerate(counterterms):
                    current = ct.nodes[0].current
                    n_subs = current.get('singular_structure').substructures
                    singular_structure = current.get('singular_structure').substructures[0]
                    all_legs = singular_structure.get_all_legs()

                    # safety check
                    if (not sorted([l.n for l in current.get('singular_structure').get_all_legs()]) == sorted([ileg,jleg,kleg])) and \
                        (not sorted([l.n for l in current.get('singular_structure').get_all_legs()]) == sorted([ileg,kleg])) and \
                        (not sorted([l.n for l in current.get('singular_structure').get_all_legs()]) == sorted([ileg,jleg])) and \
                        (not sorted([l.n for l in current.get('singular_structure').get_all_legs()]) == sorted([jleg,kleg])) and \
                        (not sorted([l.n for l in current.get('singular_structure').get_all_legs()]) == sorted([ileg])) and \
                        (not sorted([l.n for l in current.get('singular_structure').get_all_legs()]) == sorted([jleg])) and \
                        (not sorted([l.n for l in current.get('singular_structure').get_all_legs()]) == sorted([kleg])) :
                        continue

                    ### Identify cts for L1_ij
                    # L1_ij  : 3  -> [Si, Cij, SiCij]
                    # Ref. eq. 3.15, 1st line

                    if len(n_subs) == 1 and len(all_legs) == 1:
                        # Si
                        if s['sector'].id[0] == 21 :
                            s['counterterms'].append(i_ct)
                            necessary_3p_ct1_list[0] = 'S_g'
                            necessary_3p_ct1[0] = ct

                    if singular_structure.name()=='C' and len(all_legs)==2:
                        if not singular_structure.substructures:
                            # Cij
                            if sorted([l.n for l in all_legs]) == (sorted([ileg,jleg])):
                                s['counterterms'].append(i_ct)
                                if s['sector'].id[0] == 21 and s['sector'].id[1] == 21:
                                    necessary_3p_ct1_list[1] = 'C_gg'
                                elif s['sector'].id[0] == 21 and s['sector'].id[1] != 21:
                                    necessary_3p_ct1_list[1] = 'C_gq'
                                else :
                                    necessary_3p_ct1_list[1] = 'C_qqx'
                                necessary_3p_ct1[1] = ct

                    # SiCij
                    necessary_3p_ct1_list[2] = (''.join(('S_',necessary_3p_ct1_list[1])) \
                                                 if (necessary_3p_ct1_list[0] != 0 and necessary_3p_ct1_list[1] != 0) else 0)
                    necessary_3p_ct1[2] = necessary_3p_ct1[1]  # GB: why???

                    #### Identify cts for K2 ####

                    if (label == 'ijjk'):
                        # K2
                        # L2_ijjk : 7 -> [Sij,
                        #                 SCijk, SijSCijk,
                        #                 Cijk, SijCijk, SCijkCijk, SijSCijkCijk]
                        # # Ref. eq. 3.15, 2nd line

                        if len(n_subs) == 1 and len(singular_structure.substructures) == 0 :

                            if singular_structure.name()=='S' and len(all_legs)==2:
                                if not singular_structure.substructures:
                                    # Sij
                                    if sorted([l.n for l in all_legs]) == sorted([ileg,jleg]):
                                        s['counterterms'].append(i_ct)
                                        if (i_id == 21 and j_id == 21):
                                            necessary_3p_ct2_list[0] =  'SS_gg'
                                        else:
                                            necessary_3p_ct2_list[0] =  'SS_qqx'
                                        necessary_3p_ct2[0] = ct

                            if singular_structure.name()=='C' and len(all_legs)==3:
                                # Cijk
                                if sorted([l.n for l in all_legs]) == sorted([ileg,jleg,kleg]):
                                    s['counterterms'].append(i_ct)
                                    if (i_id == 21 and j_id == 21 and k_id == 21):
                                        necessary_3p_ct2_list[3] = 'CC_ggg'
                                    elif (i_id == 21 and j_id == 21 and k_id != 21):
                                        necessary_3p_ct2_list[3] = 'CC_ggq'
                                    elif (i_id == 21 and j_id != 21 and  j_id == (- k_id)):
                                        necessary_3p_ct2_list[3] = 'CC_gqqx'
                                    else:
                                        if abs(i_id) == abs(j_id) and abs(j_id) == abs(k_id):
                                            necessary_3p_ct2_list[3] = 'CC_qxqq'
                                        else:
                                            necessary_3p_ct2_list[3] = 'CC_qxqqp'
                                    necessary_3p_ct2[3] = ct

                        if len(n_subs) == 2 :
                            # here singular_structure = coll_sub -> C(i,j)
                            all_legs_C = current.get('singular_structure').substructures[1].get_all_legs()

                            # SCijk
                            if (i_id == 21 and sorted([l.n for l in all_legs_C]) == sorted([jleg,kleg])) :
                                s['counterterms'].append(i_ct)
                                if (j_id == 21 and k_id == 21):
                                    necessary_3p_ct2_list[1] = 'SC_ggg'
                                elif (j_id == 21 and k_id != 21):
                                    necessary_3p_ct2_list[1] = 'SC_ggq'
                                else:
                                    necessary_3p_ct2_list[1] = 'SC_gqqx'
                                necessary_3p_ct2[1] = ct

                        # SijSCijk
                        necessary_3p_ct2_list[2] = (''.join((necessary_3p_ct2_list[0],'_',necessary_3p_ct2_list[1])) \
                                                 if (necessary_3p_ct2_list[0] != 0 and necessary_3p_ct2_list[1] != 0) else 0)
                        necessary_3p_ct2[2] = necessary_3p_ct2[1]
                        # SijCijk
                        necessary_3p_ct2_list[4] = (''.join((necessary_3p_ct2_list[0],'_',necessary_3p_ct2_list[3])) \
                                                 if (necessary_3p_ct2_list[0] != 0 and necessary_3p_ct2_list[3] != 0) else 0)
                        necessary_3p_ct2[4] = necessary_3p_ct2[3]
                        # SCijkCijk
                        necessary_3p_ct2_list[5] = (''.join((necessary_3p_ct2_list[1],'_',necessary_3p_ct2_list[3])) \
                                                 if (necessary_3p_ct2_list[1] != 0 and necessary_3p_ct2_list[3] != 0) else 0)
                        necessary_3p_ct2[5] = necessary_3p_ct2[3]
                        # SijSCijkCijk
                        necessary_3p_ct2_list[6] = (''.join((necessary_3p_ct2_list[0],'_',necessary_3p_ct2_list[1],'_',necessary_3p_ct2_list[3])) \
                                                 if (necessary_3p_ct2_list[0] != 0 and necessary_3p_ct2_list[1] != 0 and necessary_3p_ct2_list[3] != 0) else 0)
                        necessary_3p_ct2[6] = necessary_3p_ct2[3]


                    elif (label == 'ijkj'):
                        # K2
                        # L2_ijkj : 11 -> [Sik,
                        #                  SCijk, SCkij, SikSCijk, SikSCkij,
                        #                  Cijk, SikCijk, SCijkCijk, SikSCijkCijk, SCkijCijk, SikSCkijCijk]
                        # Ref. eq. 3.15, 3rd line

                        if len(n_subs) == 1 and len(singular_structure.substructures) == 0 :

                            if singular_structure.name()=='S' and len(all_legs)==2:
                                if not singular_structure.substructures:
                                    # Sij
                                    if sorted([l.n for l in all_legs]) == sorted([ileg,kleg]):
                                        s['counterterms'].append(i_ct)
                                        if (i_id == 21 and k_id == 21):
                                            necessary_3p_ct2_list[0] =  'SS_gg'
                                        else:
                                            necessary_3p_ct2_list[0] =  'SS_qqx'
                                        necessary_3p_ct2[0] = ct

                            if singular_structure.name()=='C' and len(all_legs)==3:
                                # Cijk
                                if sorted([l.n for l in all_legs]) == sorted([ileg,jleg,kleg]):
                                    s['counterterms'].append(i_ct)
                                    if (i_id == 21 and j_id == 21 and k_id == 21):
                                        necessary_3p_ct2_list[5] = 'CC_ggg'
                                    elif (i_id == 21 and j_id == 21 and k_id != 21):
                                        necessary_3p_ct2_list[5] = 'CC_ggq'
                                    elif (i_id == 21 and j_id != 21 and  j_id == (- k_id)):
                                        necessary_3p_ct2_list[5] = 'CC_gqqx'
                                    else:
                                        if abs(i_id) == abs(j_id) and abs(j_id) == abs(k_id):
                                            necessary_3p_ct2_list[5] = 'CC_qxqq'
                                        else:
                                            necessary_3p_ct2_list[5] = 'CC_qxqqp'
                                    necessary_3p_ct2[5] = ct

                        if len(n_subs) == 2 :
                            # here singular_structure = coll_sub -> C(i,j)
                            all_legs_C = current.get('singular_structure').substructures[1].get_all_legs()

                            # SCijk
                            if (i_id == 21 and sorted([l.n for l in all_legs_C]) == sorted([jleg,kleg])) :
                                s['counterterms'].append(i_ct)
                                if (j_id == 21 and k_id == 21):
                                    necessary_3p_ct2_list[1] = 'SC_ggg'
                                elif (j_id == 21 and k_id != 21):
                                    necessary_3p_ct2_list[1] = 'SC_ggq'
                                else:
                                    necessary_3p_ct2_list[1] = 'SC_gqqx'
                                necessary_3p_ct2[1] = ct
                            # SCkij
                            if (k_id == 21 and sorted([l.n for l in all_legs_C]) == sorted([ileg,jleg])) :
                                s['counterterms'].append(i_ct)
                                if (i_id == 21 and j_id == 21):
                                    necessary_3p_ct2_list[2] = 'SC_ggg'
                                elif (i_id == 21 and j_id != 21):
                                    necessary_3p_ct2_list[2] = 'SC_ggq'
                                else:
                                    necessary_3p_ct2_list[2] = 'SC_gqqx'
                                necessary_3p_ct2[2] = ct

                        # SikSCijk
                        necessary_3p_ct2_list[3] = (''.join((necessary_3p_ct2_list[0],'_',necessary_3p_ct2_list[1])) \
                                                 if (necessary_3p_ct2_list[0] != 0 and necessary_3p_ct2_list[1] != 0) else 0)
                        necessary_3p_ct2[3] = necessary_3p_ct2[1]
                        # SikSCkij
                        necessary_3p_ct2_list[4] = (''.join((necessary_3p_ct2_list[0],'_',necessary_3p_ct2_list[2])) \
                                                 if (necessary_3p_ct2_list[0] != 0 and necessary_3p_ct2_list[2] != 0) else 0)
                        necessary_3p_ct2[4] = necessary_3p_ct2[2]
                        # SikCijk
                        necessary_3p_ct2_list[6] = (''.join((necessary_3p_ct2_list[0],'_',necessary_3p_ct2_list[5])) \
                                                 if (necessary_3p_ct2_list[0] != 0 and necessary_3p_ct2_list[5] != 0) else 0)
                        necessary_3p_ct2[6] = necessary_3p_ct2[5]
                        # SCijkCijk
                        necessary_3p_ct2_list[7] = (''.join((necessary_3p_ct2_list[1],'_',necessary_3p_ct2_list[5])) \
                                                 if (necessary_3p_ct2_list[1] != 0 and necessary_3p_ct2_list[5] != 0) else 0)
                        necessary_3p_ct2[7] = necessary_3p_ct2[5]
                        # SikSCijkCijk
                        necessary_3p_ct2_list[8] = (''.join((necessary_3p_ct2_list[0],'_',necessary_3p_ct2_list[1],'_',necessary_3p_ct2_list[5])) \
                                                 if (necessary_3p_ct2_list[0] != 0 and necessary_3p_ct2_list[1] != 0 and necessary_3p_ct2_list[5] != 0) else 0)
                        necessary_3p_ct2[8] = necessary_3p_ct2[5]
                        # SCkijCijk
                        necessary_3p_ct2_list[9] = (''.join((necessary_3p_ct2_list[2],'_',necessary_3p_ct2_list[5])) \
                                                 if (necessary_3p_ct2_list[2] != 0 and necessary_3p_ct2_list[5] != 0) else 0)
                        necessary_3p_ct2[9] = necessary_3p_ct2[5]
                        # SikSCkijCijk
                        necessary_3p_ct2_list[10] = (''.join((necessary_3p_ct2_list[0],'_',necessary_3p_ct2_list[2],'_',necessary_3p_ct2_list[5])) \
                                                 if (necessary_3p_ct2_list[0] != 0 and necessary_3p_ct2_list[2] != 0 and necessary_3p_ct2_list[5] != 0) else 0)
                        necessary_3p_ct2[10] = necessary_3p_ct2[5]


                    #### Identify cts for K12 ####

                    if (label == 'ijjk'):
                        # K12
                        # L12_ijjk : 13 -> [Si Sij, Si SCijk, Si SijSCijk,
                        #                   Si Cijk, Si SijCijk, Si SCijkCijk, Si SijSCijkCijk,
                        #                   Cij Sij, Cij Cijk, Cij SijCijk,
                        #                   SiCij Sij, SiCij Cijk, SiCij SijCijk]
                        # Ref. eq. 3.15, 5th line

                        # Si Sij
                        necessary_3p_ct12_list[0] = (''.join(('S_',necessary_3p_ct2_list[0])) \
                                                 if (necessary_3p_ct1_list[0] != 0 and necessary_3p_ct2_list[0] != 0) else 0)
                        # Si SCijk
                        necessary_3p_ct12_list[1] = (''.join(('S_',necessary_3p_ct2_list[1])) \
                                                 if (necessary_3p_ct1_list[0] != 0 and necessary_3p_ct2_list[1] != 0) else 0)
                        # Si SijSCijk
                        necessary_3p_ct12_list[2] = (''.join(('S_',necessary_3p_ct2_list[2])) \
                                                 if (necessary_3p_ct1_list[0] != 0 and necessary_3p_ct2_list[2] != 0) else 0)
                        # Si Cijk
                        necessary_3p_ct12_list[3] = (''.join(('S_',necessary_3p_ct2_list[3])) \
                                                 if (necessary_3p_ct1_list[0] != 0 and necessary_3p_ct2_list[3] != 0) else 0)
                        # Si SijCijk
                        necessary_3p_ct12_list[4] = (''.join(('S_',necessary_3p_ct2_list[4])) \
                                                 if (necessary_3p_ct1_list[0] != 0 and necessary_3p_ct2_list[4] != 0) else 0)
                        # Si SCijkCijk
                        necessary_3p_ct12_list[5] = (''.join(('S_',necessary_3p_ct2_list[5])) \
                                                 if (necessary_3p_ct1_list[0] != 0 and necessary_3p_ct2_list[5] != 0) else 0)
                        # Si SijSCijkCijk
                        necessary_3p_ct12_list[6] = (''.join(('S_',necessary_3p_ct2_list[6])) \
                                                 if (necessary_3p_ct1_list[0] != 0 and necessary_3p_ct2_list[6] != 0) else 0)
                        # Cij Sij
                        necessary_3p_ct12_list[7] = (''.join(('C_',necessary_3p_ct2_list[0])) \
                                                 if (necessary_3p_ct1_list[1] != 0 and necessary_3p_ct2_list[0] != 0) else 0)
                        # Cij Cijk
                        necessary_3p_ct12_list[8] = (''.join(('C_',necessary_3p_ct2_list[3])) \
                                                 if (necessary_3p_ct1_list[1] != 0 and necessary_3p_ct2_list[3] != 0) else 0)
                        # Cij SijCijk
                        necessary_3p_ct12_list[9] = (''.join(('C_',necessary_3p_ct2_list[4])) \
                                                 if (necessary_3p_ct1_list[1] != 0 and necessary_3p_ct2_list[4] != 0) else 0)
                        # SiCij Sij
                        necessary_3p_ct12_list[10] = (''.join(('SC_',necessary_3p_ct2_list[0])) \
                                                 if (necessary_3p_ct1_list[0] != 0 and necessary_3p_ct1_list[1] != 0 and necessary_3p_ct2_list[0] != 0) else 0)
                        # SiCij Cijk
                        necessary_3p_ct12_list[11] = (''.join(('SC_',necessary_3p_ct2_list[3])) \
                                                 if (necessary_3p_ct1_list[0] != 0 and necessary_3p_ct1_list[1] != 0 and necessary_3p_ct2_list[3] != 0) else 0)
                        # SiCij SijCijk
                        necessary_3p_ct12_list[12] = (''.join(('SC_',necessary_3p_ct2_list[4])) \
                                                 if (necessary_3p_ct1_list[0] != 0 and necessary_3p_ct1_list[1] != 0 and necessary_3p_ct2_list[4] != 0) else 0)

                    elif (label == 'ijkj'):
                        # K12
                        # L12_ijkj : 13 -> [Si Sik, Si SCijk, Si SikSCijk,
                        #                   Si Cijk, Si SikCijk, Si SCijkCijk, Si SikSCijkCijk,
                        #                   Cij SCkij, Cij Cijk, Cij SCkijCijk,
                        #                   SiCij SCkij, SiCij Cijk, SiCij SCkijCijk]
                        # Ref. eq. 3.15, 6th line

                        # Si Sik
                        necessary_3p_ct12_list[0] = (''.join(('S_',necessary_3p_ct2_list[0])) \
                                                 if (necessary_3p_ct1_list[0] != 0 and necessary_3p_ct2_list[0] != 0) else 0)
                        # Si SCijk
                        necessary_3p_ct12_list[1] = (''.join(('S_',necessary_3p_ct2_list[1])) \
                                                 if (necessary_3p_ct1_list[0] != 0 and necessary_3p_ct2_list[1] != 0) else 0)
                        # Si SikSCijk
                        necessary_3p_ct12_list[2] = (''.join(('S_',necessary_3p_ct2_list[3])) \
                                                 if (necessary_3p_ct1_list[0] != 0 and necessary_3p_ct2_list[3] != 0) else 0)
                        # Si Cijk
                        necessary_3p_ct12_list[3] = (''.join(('S_',necessary_3p_ct2_list[5])) \
                                                 if (necessary_3p_ct1_list[0] != 0 and necessary_3p_ct2_list[5] != 0) else 0)
                        # Si SikCijk
                        necessary_3p_ct12_list[4] = (''.join(('S_',necessary_3p_ct2_list[6])) \
                                                 if (necessary_3p_ct1_list[0] != 0 and necessary_3p_ct2_list[6] != 0) else 0)
                        # Si SCijkCijk
                        necessary_3p_ct12_list[5] = (''.join(('S_',necessary_3p_ct2_list[7])) \
                                                 if (necessary_3p_ct1_list[0] != 0 and necessary_3p_ct2_list[7] != 0) else 0)
                        # Si SikSCijkCijk
                        necessary_3p_ct12_list[6] = (''.join(('S_',necessary_3p_ct2_list[8])) \
                                                 if (necessary_3p_ct1_list[0] != 0 and necessary_3p_ct2_list[8] != 0) else 0)
                        # Cij SCkij
                        necessary_3p_ct12_list[7] = (''.join(('C_',necessary_3p_ct2_list[2])) \
                                                 if (necessary_3p_ct1_list[1] != 0 and necessary_3p_ct2_list[2] != 0) else 0)
                        # Cij Cijk
                        necessary_3p_ct12_list[8] = (''.join(('C_',necessary_3p_ct2_list[5])) \
                                                 if (necessary_3p_ct1_list[1] != 0 and necessary_3p_ct2_list[5] != 0) else 0)
                        # Cij SCkijCijk
                        necessary_3p_ct12_list[9] = (''.join(('C_',necessary_3p_ct2_list[9])) \
                                                 if (necessary_3p_ct1_list[1] != 0 and necessary_3p_ct2_list[9] != 0) else 0)
                        # SiCij SCkij
                        necessary_3p_ct12_list[10] = (''.join(('SC_',necessary_3p_ct2_list[2])) \
                                                 if (necessary_3p_ct1_list[0] != 0 and necessary_3p_ct1_list[1] != 0 and necessary_3p_ct2_list[2] != 0) else 0)
                        # SiCij Cijk
                        necessary_3p_ct12_list[11] = (''.join(('SC_',necessary_3p_ct2_list[5])) \
                                                 if (necessary_3p_ct1_list[0] != 0 and necessary_3p_ct1_list[1] != 0 and necessary_3p_ct2_list[5] != 0) else 0)
                        # SiCij SCkijCijk
                        necessary_3p_ct12_list[12] = (''.join(('SC_',necessary_3p_ct2_list[9])) \
                                                 if (necessary_3p_ct1_list[0] != 0 and necessary_3p_ct1_list[1] != 0 and necessary_3p_ct2_list[9] != 0) else 0)


                # Collect all sector information
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


        ######## 4p #########
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
                necessary_4p_ct1_list = [0] * (3)
                necessary_4p_ct1 = [0] * (3)
                necessary_4p_ct2_list = [0] * (9)
                necessary_4p_ct2 = [0] * (9)
                necessary_4p_ct12_list = [0] * (9)

                print('****** NEW SECTOR ******')
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
                    i_id = s['sector'].id[0]
                    j_id = s['sector'].id[1]
                    k_id = s['sector'].id[2]
                    l_id = s['sector'].id[3]

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


                    # L1_ijkl  : 3 -> [Si, Cij, SiCij]

                    if len(n_subs) == 1 and len(all_legs) == 1:
                        # Si
                        if i_id == 21 :
                            s['counterterms'].append(i_ct)
                            necessary_4p_ct1_list[0] = 'S_g'
                            necessary_4p_ct1[0] = ct

                    if singular_structure.name()=='C' and len(all_legs)==2:
                        if not singular_structure.substructures:
                            # Cij
                            if sorted([l.n for l in all_legs]) == (sorted([ileg,jleg])):
                                s['counterterms'].append(i_ct)
                                if (i_id == 21 and j_id == 21):
                                    necessary_4p_ct1_list[1] = 'C_gg'
                                elif (i_id == 21 and j_id != 21):
                                    necessary_4p_ct1_list[1] = 'C_gq'
                                else :
                                    necessary_4p_ct1_list[1] = 'C_qqx'
                                necessary_4p_ct1[1] = ct
                            # Ckl
                            Ckl_flag = False
                            if sorted([l.n for l in all_legs]) == (sorted([kleg,lleg])):
                                Ckl_flag = True

                    # SiCij
                    necessary_4p_ct1_list[2] = (''.join(('S_',necessary_4p_ct1_list[1])) \
                                                 if (necessary_4p_ct1_list[0] != 0 and necessary_4p_ct1_list[1] != 0) else 0)
                    necessary_4p_ct1[2] = necessary_4p_ct1[1]


                    # L2_ijkl  : 9  ->  [Sik, SCikl, SCkij, SikSCikl, SikSCkij,
                    #                    Cijkl, SikCijkl, SCiklCijkl, SCkijCijkl]

                    if len(n_subs) == 1 and len(singular_structure.substructures) == 0 :

                        if singular_structure.name()=='S' and len(all_legs)==2:
                            # Sik
                            if sorted([l.n for l in all_legs]) == sorted([ileg,kleg]):
                                s['counterterms'].append(i_ct)
                                if (i_id == 21 and k_id == 21):
                                    necessary_4p_ct2_list[0] =  'SS_gg'
                                else:
                                    necessary_4p_ct2_list[0] =  'SS_qqx'
                                necessary_4p_ct2[0] = ct

                    if len(n_subs) == 2 and len(singular_structure.substructures) == 0 :

                        if n_subs[0].name()=='S' and n_subs[1].name()=='C':
                            # SCikl
                            if i_id == 21 and sorted([l.n for l in n_subs[1].get_all_legs()]) == (sorted([kleg,lleg])):
                                if (k_id == 21 and l_id == 21):
                                    necessary_4p_ct2_list[1] = 'SC_ggg'
                                elif (k_id == 21 and l_id != 21):
                                    necessary_4p_ct2_list[1] = 'SC_ggq'
                                else:
                                    necessary_4p_ct2_list[1] = 'SC_gqqx'
                                necessary_4p_ct2[1] = ct
                            # SCkij
                            if k_id == 21 and sorted([l.n for l in n_subs[1].get_all_legs()]) == (sorted([ileg,jleg])):
                                if (i_id == 21 and j_id == 21):
                                    necessary_4p_ct2_list[2] = 'SC_ggg'
                                elif (i_id == 21 and j_id != 21):
                                    necessary_4p_ct2_list[2] = 'SC_ggq'
                                else:
                                    necessary_4p_ct2_list[2] = 'SC_gqqx'
                                necessary_4p_ct2[2] = ct

                        # Cijkl
                        if n_subs[0].name()=='C' and n_subs[1].name()=='C':
                            if sorted([l.n for l in n_subs[0].get_all_legs()]) == (sorted([ileg,jleg])) and \
                                sorted([l.n for l in n_subs[1].get_all_legs()]) == (sorted([kleg,lleg])) and \
                                    Ckl_flag == True:
                                if i_id == 21 and j_id == 21 and k_id == 21 and l_id == 21:
                                    necessary_4p_ct2_list[5] = 'CC_gggg'
                                elif i_id == 21 and j_id == 21 and k_id == 21 and l_id != 21:
                                    necessary_4p_ct2_list[5] = 'CC_gggq'
                                elif k_id == 21 and l_id == 21 and i_id == 21 and j_id != 21:
                                    necessary_4p_ct2_list[5] = 'CC_gqgg'
                                elif i_id == 21 and j_id == 21 and abs(k_id) == abs(l_id):
                                    necessary_4p_ct2_list[5] = 'CC_ggqxq'
                                elif k_id == 21 and l_id == 21 and abs(i_id) == abs(j_id):
                                    necessary_4p_ct2_list[5] = 'CC_qxqgg'
                                elif i_id == 21 and j_id != 21 and k_id == 21 and l_id != 21:
                                    necessary_4p_ct2_list[5] = 'CC_gqgq'
                                elif i_id == 21 and j_id != 21 and abs(k_id) == abs(l_id):
                                    necessary_4p_ct2_list[5] = 'CC_gqqxq'
                                elif k_id == 21 and l_id != 21 and abs(i_id) == abs(j_id):
                                    necessary_4p_ct2_list[5] = 'CC_qxqgq'
                                elif abs(i_id) == abs(j_id) and abs(k_id) == abs(l_id):
                                    necessary_4p_ct2_list[5] = 'CC_qxqqxq'
                                necessary_4p_ct2[5] = ct

                    # SikSCikl
                    necessary_4p_ct2_list[3] = (''.join((necessary_4p_ct2_list[0],'_',necessary_4p_ct2_list[1])) \
                                                 if (necessary_4p_ct2_list[0] != 0 and necessary_4p_ct2_list[1] != 0) else 0)
                    necessary_4p_ct2[3] = necessary_4p_ct2[1]
                    # SikSCkij
                    necessary_4p_ct2_list[4] = (''.join((necessary_4p_ct2_list[0],'_',necessary_4p_ct2_list[2])) \
                                                 if (necessary_4p_ct2_list[0] != 0 and necessary_4p_ct2_list[2] != 0) else 0)
                    necessary_4p_ct2[4] = necessary_4p_ct2[2]
                    # SikCijkl
                    necessary_4p_ct2_list[6] = (''.join((necessary_4p_ct2_list[0],'_',necessary_4p_ct2_list[5])) \
                                                 if (necessary_4p_ct2_list[0] != 0 and necessary_4p_ct2_list[5] != 0) else 0)
                    necessary_4p_ct2[6] = necessary_4p_ct2[5]
                    # SCiklCijkl
                    necessary_4p_ct2_list[7] = (''.join((necessary_4p_ct2_list[1],'_',necessary_4p_ct2_list[5])) \
                                                 if (necessary_4p_ct2_list[1] != 0 and necessary_4p_ct2_list[5] != 0) else 0)
                    necessary_4p_ct2[7] = necessary_4p_ct2[5]
                    # SCkijCijkl
                    necessary_4p_ct2_list[8] = (''.join((necessary_4p_ct2_list[2],'_',necessary_4p_ct2_list[5])) \
                                                 if (necessary_4p_ct2_list[2] != 0 and necessary_4p_ct2_list[5] != 0) else 0)
                    necessary_4p_ct2[8] = necessary_4p_ct2[5]

                    # L12_ijkl : 9  ->  [Si Sik, Si SCikl, Si SikSCikl
                    #                    Cij SCkij, Cij Cijkl, Cij SCkijCijkl,
                    #                    SiCij SCkij, SiCij Cijkl, SiCij SCkijCijkl]

                    # Si Sik
                    necessary_4p_ct12_list[0] = (''.join(('S_',necessary_4p_ct2_list[0])) \
                                                 if (necessary_4p_ct1_list[0] != 0 and necessary_4p_ct2_list[0] != 0) else 0)
                    # Si SCikl
                    necessary_4p_ct12_list[1] = (''.join(('S_',necessary_4p_ct2_list[1])) \
                                                 if (necessary_4p_ct1_list[0] != 0 and necessary_4p_ct2_list[1] != 0) else 0)
                    # Si SikSCikl
                    necessary_4p_ct12_list[2] = (''.join(('S_',necessary_4p_ct2_list[3])) \
                                                 if (necessary_4p_ct1_list[0] != 0 and necessary_4p_ct2_list[3] != 0) else 0)
                    # Cij SCkij
                    necessary_4p_ct12_list[3] = (''.join(('C_',necessary_4p_ct2_list[2])) \
                                                 if (necessary_4p_ct1_list[1] != 0 and necessary_4p_ct2_list[2] != 0) else 0)
                    # Cij Cijkl
                    necessary_4p_ct12_list[4] = (''.join(('C_',necessary_4p_ct2_list[5])) \
                                                 if (necessary_4p_ct1_list[1] != 0 and necessary_4p_ct2_list[5] != 0) else 0)
                    # Cij SCkijCijkl
                    necessary_4p_ct12_list[5] = (''.join(('C_',necessary_4p_ct2_list[8])) \
                                                 if (necessary_4p_ct1_list[1] != 0 and necessary_4p_ct2_list[8] != 0) else 0)
                    # SiCij SCkij
                    necessary_4p_ct12_list[6] = (''.join(('SC_',necessary_4p_ct2_list[2])) \
                                                 if (necessary_4p_ct1_list[0] != 0 and necessary_4p_ct1_list[1] != 0 and necessary_4p_ct2_list[2] != 0) else 0)
                    # SiCij Cijkl
                    necessary_4p_ct12_list[7] = (''.join(('SC_',necessary_4p_ct2_list[5])) \
                                                 if (necessary_4p_ct1_list[0] != 0 and necessary_4p_ct1_list[1] != 0 and necessary_4p_ct2_list[5] != 0) else 0)
                    # SiCij SCkijCijkl
                    necessary_4p_ct12_list[8] = (''.join(('SC_',necessary_4p_ct2_list[8])) \
                                                 if (necessary_4p_ct1_list[0] != 0 and necessary_4p_ct1_list[1] != 0 and necessary_4p_ct2_list[8] != 0) else 0)

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
        dirpath = pjoin(dirpath, 'SubProcesses', "P%s" % defining_process.shell_string())
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

# ######### Write NNLO_K_isec_jsec_jsec_ksec.f and NNLO_K_isec_jsec_ksec_jsec.f (3-particle sector)

        # Set replace_dict for NNLO_K
        replace_dict_ct = {}
        replace_dict_limits = {}
        replace_dict_double_real = {}
        necessary_default_3p_ct_list = ['S_g', 'C_gg', 'C_gq', 'C_qqx', 'S_C_gg', 'S_C_gq', \
                                        'SS_gg', 'SS_qqx', \
                                        'SC_ggg', 'SC_ggq', 'SC_gqqx', \
                                        'SS_gg_SC_ggg', 'SS_gg_SC_ggq', \
                                        'CC_ggg', 'CC_ggq', 'CC_gqqx', 'CC_qxqq', 'CC_qxqqp', \
                                        'SS_gg_CC_ggg', 'SS_gg_CC_ggq', 'SS_qqx_CC_qxqq', 'SS_qqx_CC_qxqqp', \
                                        'SC_ggg_C_ggg', 'SC_ggq_CC_ggq', 'SC_gqqx_CC_gqqx', \
                                        'SS_gg_SC_ggg_C_ggg', 'SS_gg_SC_ggq_CC_ggq', \
                                        'S_SS_gg', 'S_SC_ggg', 'S_SC_ggq', 'S_SC_gqqx', \
                                        'S_SS_gg_SC_ggg', 'S_SS_gg_SC_ggq', \
                                        'S_CC_ggg', 'S_CC_ggq', 'S_CC_gqqx', \
                                        'S_SS_gg_CC_ggg', 'S_SS_gg_CC_ggq', \
                                        'S_SC_ggg_C_ggg', 'S_SC_ggq_CC_ggq', 'S_SC_gqqx_CC_gqqx', \
                                        'S_SS_gg_SC_ggg_C_ggg', 'S_SS_gg_SC_ggq_CC_ggq', \
                                        'C_SS_gg', 'C_SS_qqx', \
                                        'C_SC_ggg', 'C_SC_ggq', 'C_SC_gqqx', \
                                        'C_CC_ggg', 'C_CC_ggq', 'C_CC_gqqx', 'C_CC_qxqq', 'C_CC_qxqqp', \
                                        'C_SS_gg_CC_ggg', 'C_SS_gg_CC_ggq', 'C_SS_qqx_CC_qxqq', 'C_SS_qqx_CC_qxqqp', \
                                        'C_SC_ggg_C_ggg', 'C_SC_ggq_CC_ggq', 'C_SC_gqqx_CC_gqqx', \
                                        'SC_SS_gg', \
                                        'SC_SC_ggg', 'SC_SC_ggq', 'SC_SC_gqqx', \
                                        'SC_CC_ggg', 'SC_CC_ggq', 'SC_CC_gqqx', \
                                        'SC_SS_gg_CC_ggg', 'SC_SS_gg_CC_ggq', \
                                        'SC_SC_ggg_C_ggg', 'SC_SC_ggq_CC_ggq', 'SC_SC_gqqx_CC_gqqx']

        # K2_ijjk  : 7 -> [Sij, SCijk, SijSCijk, Cijk, SijCijk, SCijkCijk, SijSCijkCijk]
        # K2_ijkj  : 11 -> [Sik, SCijk, SCkij, SikSCijk, SikSCkij, Cijk, SikCijk, SCijkCijk, SikSCijkCijk, SCkijCijk, SikSCkijCijk]
        # K12_ijjk : 13 -> [Si Sij, Si SCijk, Si SijSCijk, Si Cijk, Si SijCijk, Si SCijkCijk, Si SijSCijkCijk,
        #                   Cij Sij, Cij Cijk, Cij SijCijk, SiCij Sij, SiCij Cijk, SiCij SijCijk]
        # K12_ijkj : 13 -> [Si Sik, Si SCijk, Si SikSCijk, Si Cijk, Si SikCijk, Si SCijkCijk, Si SikSCijkCijk,
        #                   Cij SCkij, Cij Cijk, Cij SCkijCijk, SiCij SCkij, SiCij Cijk, SiCij SCkijCijk]
        K1_labels = ['S_i', 'C_ij', 'S_i C_ij']
        K2_labels_ijjk = ['SS_ij', 'SC_ijk', 'SS_ij SC_ijk', 'CC_ijk', 'SS_ij CC_ijk', 'SC_ijk CC_ijk', 'SS_ij SC_ijk CC_ijk']
        K2_labels_ijkj = ['SS_ik', 'SC_ijk', 'SC_kij', 'SS_ik SC_ijk', 'SS_ik SCkij', \
                          'CC_ijk', 'SS_ik CC_ijk', 'SC_ijk CC_ijk', 'SS_ik SC_ijk CC_ijk', 'SC_kij CC_ijk', 'SS_ik SC_kij CC_ijk']
        K12_labels_ijjk = ['S_i SS_ij', 'S_i SC_ijk', 'S_i SS_ij SC_ijk', \
                           'S_i CC_ijk', 'S_i SS_ij CC_ijk', 'S_i SC_ijk CC_ijk', 'S_i SS_ij SC_ijk CC_ijk', \
                           'C_ij SS_ij', 'C_ij CC_ijk', 'C_ij SS_ij CC_ijk', \
                           'S_i C_ij SS_ij', 'S_i C_ij CC_ijk', 'S_i C_ij SS_ij CC_ijk']
        K12_labels_ijkj = ['S_i SS_ik', 'S_i SC_ijk', 'S_i SS_ik SC_ijk', \
                           'S_i CC_ijk', 'S_i SS_ik CC_ijk', 'S_i SC_ijk CC_ijk', 'S_i SS_ik SC_ijk CC_ijk', \
                           'C_ij SC_kij', 'C_ij CC_ijk', 'C_ij SC_kij CC_ijk', \
                           'S_i C_ij SC_kij', 'S_i C_ij CC_ijk', 'S_i C_ij SC_kij CC_ijk']

        # Rule: examples
        # Cijk      -> isec,jsec,ksec
        # Sik SCijk -> isec,ksec,jsec
        # Cij SCkij -> ksec,isec,jsec
        K1_3p_indices = ['isec', 'isec,jsec', 'isec,jsec']
        K2_3p_indices_ijjk = ['isec,jsec', 'isec,jsec,ksec', 'isec,jsec,ksec', 'isec,jsec,ksec', 'isec,jsec,ksec', 'isec,jsec,ksec', 'isec,jsec,ksec']
        K2_3p_indices_ijkj = ['isec,ksec', 'isec,jsec,ksec', 'ksec,isec,jsec', 'isec,ksec,jsec', 'ksec,isec,jsec', \
                              'isec,jsec,ksec', 'isec,ksec,jsec', 'isec,jsec,ksec', 'isec,ksec,jsec', 'ksec,isec,jsec', 'ksec,isec,jsec']
        K12_3p_indices_ijjk = ['isec,jsec', 'isec,jsec,ksec', 'isec,jsec,ksec', 'isec,jsec,ksec', 'isec,jsec,ksec', 'isec,jsec,ksec', 'isec,jsec,ksec', \
                               'isec,jsec', 'isec,jsec,ksec', 'isec,jsec,ksec', 'isec,jsec', 'isec,jsec,ksec', 'isec,jsec,ksec']
        K12_3p_indices_ijkj = ['isec,ksec', 'isec,jsec,ksec', 'isec,ksec,jsec', 'isec,jsec,ksec', 'isec,ksec,jsec', 'isec,jsec,ksec', 'isec,ksec,jsec', \
                               'ksec,isec,jsec', 'isec,jsec,ksec', 'ksec,isec,jsec', \
                               'ksec,isec,jsec', 'isec,jsec,ksec', 'ksec,isec,jsec']

        label = ''
        for i in range(0,len(all_3p_sector_list)):
            list_str_defK1 = []
            list_str_M2_K1 = []
            list_str_defK2 = []
            list_str_M2_K2 = []
            list_str_defK12 = []
            list_str_M2_K12 = []
            mapping = []
            lsec = 0    # no fourth index needed
            if (all_3p_sector_list[i][1] == all_3p_sector_list[i][2]):
                label = 'ijjk'
                isec = all_3p_sector_list[i][0]
                jsec = all_3p_sector_list[i][1]
                ksec = all_3p_sector_list[i][3]
                id_isec = all_3p_sector_id_list[i][0]
                id_jsec = all_3p_sector_id_list[i][1]
                id_ksec = all_3p_sector_id_list[i][3]
                K2_labels = K2_labels_ijjk
                K12_labels = K12_labels_ijjk
                K2_3p_indices = K2_3p_indices_ijjk
                K12_3p_indices = K12_3p_indices_ijjk
            elif (all_3p_sector_list[i][1] == all_3p_sector_list[i][3]):
                label = 'ijkj'
                isec = all_3p_sector_list[i][0]
                jsec = all_3p_sector_list[i][1]
                ksec = all_3p_sector_list[i][2]
                id_isec = all_3p_sector_id_list[i][0]
                id_jsec = all_3p_sector_id_list[i][1]
                id_ksec = all_3p_sector_id_list[i][2]
                K2_labels = K2_labels_ijkj
                K12_labels = K12_labels_ijkj
                K2_3p_indices = K2_3p_indices_ijkj
                K12_3p_indices = K12_3p_indices_ijkj

            # Extract the reference particle leg from recoiler_function.py
            iref = all_3p_sector_recoilers[i]
            replace_dict_ct['iref'] = iref
            if (isec == iref) or (jsec == iref) or (ksec == iref):
                        raise MadEvent7Error('Wrong recoiler %d,%d,%d,%d!' % (isec,jsec,ksec,iref))
            # Check sector indicies
            if isec == jsec or isec == ksec or jsec == ksec:
                raise MadEvent7Error('Wrong sector indices %d,%d,%d!' % (isec,jsec,ksec))

            replace_dict_ct['isec'] = isec
            replace_dict_ct['jsec'] = jsec
            replace_dict_ct['ksec'] = ksec
            replace_dict_ct['lsec'] = lsec              # trivial = 0
            replace_dict_double_real['isec'] = isec
            replace_dict_double_real['jsec'] = jsec
            replace_dict_double_real['ksec'] = ksec
            replace_dict_double_real['iref'] = iref
            replace_dict_double_real['lsec'] = lsec     # trivial = 0
            replace_dict_limits['isec'] = isec
            replace_dict_limits['jsec'] = jsec
            replace_dict_limits['ksec'] = ksec
            replace_dict_limits['lsec'] = lsec          # trivial = 0
            # Write sector indices in NNLO_RRsub
            if label == 'ijjk':
                replace_dict_double_real['c3p'] = jsec
                replace_dict_double_real['d3p'] = ksec
                replace_dict_ct['c3p'] = jsec
                replace_dict_ct['d3p'] = ksec
            elif label == 'ijkj':
                replace_dict_double_real['c3p'] = ksec
                replace_dict_double_real['d3p'] = jsec
                replace_dict_ct['c3p'] = ksec
                replace_dict_ct['d3p'] = jsec

            replace_dict_limits['proc_prefix_rr'] = str(defining_process.shell_string(schannel=True,
                                        forbid=True, main=False, pdg_order=False, print_id = False))
            replace_dict_double_real['proc_prefix_rr'] = str(defining_process.shell_string(schannel=True,
                                       forbid=True, main=False, pdg_order=False, print_id = False) + '_')
            replace_dict_double_real['str_UBorn'] = 'dummy'
            replace_dict_double_real['UBgraphs'] = 'dummy'

            # Update sector_info dictionary
            sector_info = {
                'isec'          :   0,
                'jsec'          :   0,
                'ksec'          :   0,
                'lsec'          :   0,
                'iref'          :   0,
                'c3p'           :   0,
                'd3p'           :   0,
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
            if label == 'ijjk':
                sector_info['c3p'] = jsec
                sector_info['d3p'] = ksec
            elif label == 'ijkj':
                sector_info['c3p'] = ksec
                sector_info['d3p'] = jsec

            # default mapping
            mapping = [('isec', isec), ('jsec', jsec), ('ksec', ksec), ('iref', iref)]
            sector_info['mapping'] = [mapping[0][1], mapping[1][1], mapping[2][1], mapping[3][1]]
            # specify ((isec,jsec,iref),(jsec,ksec,iref)) mapping choice
            mapping_str = """ \
                iU1 = %s
                iS1 = %s
                iB1 = %s
                iU2 = %s
                iS2 = %s
                iB2 = %s
                iA1 = 1 ! default azimuth for NLO
                iA2 = 1 ! default azimuth for NLO
            """ % (mapping[0][0], mapping[1][0], mapping[3][0], mapping[1][0],mapping[2][0],mapping[3][0])


            # Initialise NLO_IR_limits.f for every sector [ij]
            list = []
            if label == 'ijjk':
                list.append('c Collection of relevant limits for sector [%d,%d,%d,%d] \n' %(isec,jsec,jsec,ksec))
            elif label == 'ijkj':
                list.append('c Collection of relevant limits for sector [%d,%d,%d,%d] \n' %(isec,jsec,ksec,jsec))
            list.append('c K1  ' + str(all_3p_K1_ct[i])  + ' \n')
            list.append('c K2  ' + str(all_3p_K2_ct[i])  + ' \n')
            list.append('c K12 ' + str(all_3p_K12_ct[i]) + ' \n')
            string = "".join(list)
            NNLO_IR_limits_tmp_path = dirmadnklo + '/tmp_fortran/tmp_files/NNLO_limits/'
            with open(NNLO_IR_limits_tmp_path + "IR_tmp.f", "w") as f:
                f.writelines(string)

            # Loop on K1 cts
            # Recall that all_3p_K1_ct = [Si, Cij, SiCij]
            ct_list = []
            write_S = True
            write_HC = True
            for j in range(0, len(all_3p_K1_ct[i])):
                if all_3p_K1_ct[i][j] ==  0:
                    continue
                if j == 0: # S_i
                    list_str_M2_K1.append('K%s=K%s+M2_%s(%s,xs,xp,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                       % (all_3p_K1_ct[i][j].split("_")[0], all_3p_K1_ct[i][j].split("_")[0], all_3p_K1_ct[i][j], K1_3p_indices[j]))
                    list_str_M2_K1.append('if(ierr.eq.1)goto 999\n')
                    # This following if statement is probabibly no longer necessary as in the W sector the sector can anyhow appear once in the collection of tmp kernels
                    if (write_S): # just one call to the kernel is needed
                        write_S = False
                        os.system('cat ' + NNLO_IR_limits_tmp_path + all_3p_K1_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
                elif j == 1: # C_ij
                    list_str_M2_K1.append('K%s=K%s+M2_%s(%s,iref,xs,xp,xsb,xpb,wgt,xj,nitRR,1d0,wgt_chan,ierr)\n'
                                       % (all_3p_K1_ct[i][j].split("_")[0], all_3p_K1_ct[i][j].split("_")[0], all_3p_K1_ct[i][j], K1_3p_indices[j]))
                    list_str_M2_K1.append('if(ierr.eq.1)goto 999\n')
                    if (write_HC): # just one call to the kernel is needed
                        write_HC = False
                        os.system('cat ' + NNLO_IR_limits_tmp_path + all_3p_K1_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
                elif j == 2: # S_C_ij
                    list_str_M2_K1.append('K%s=K%s+M2_%s(%s,xs,xp,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                       % ('SC', 'SC', all_3p_K1_ct[i][j], K1_3p_indices[j]))
                    list_str_M2_K1.append('if(ierr.eq.1)goto 999\n')
                    os.system('cat ' + NNLO_IR_limits_tmp_path + all_3p_K1_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')

                # Extract underlying real string
                self.get_uproc_str('Real', uB_all_3p_K1_ct[i][j], all_3p_K1_ct[i][j], dirpathR_head, replace_dict_limits,
                                       replace_dict_double_real, UReal_procs, path_UReal_procs, sector_info)

                if all_3p_K1_ct[i][j] not in ct_list:
                    ct_list.append(all_3p_K1_ct[i][j])
                    tmp_str = """
c       %s
        DOUBLE PRECISION M2_%s""" %(K1_labels[j],all_3p_K1_ct[i][j])
                    list_str_defK1.append(tmp_str)


            # Loop on K2 cts
            # L2_ijjk : 7 -> [Sij, SCijk, SijSCijk, Cijk, SijCijk, SCijkCijk, SijSCijkCijk]

            ct_list = []
            if label == 'ijjk':
                tmp_str = """
c       KSS = SS_ij
c       KSC = SC_ijk (1 - SS_ij)
c       KCC = CC_ijk (1 - SS_ij) (1 - SC_ijk)"""
                list_str_defK2.append(tmp_str)
                for j in range(0, len(all_3p_K2_ct[i])):
                    if all_3p_K2_ct[i][j] ==  0:
                        continue
                    if j == 0: # SS_ij
                        list_str_M2_K2.append('K%s=K%s+M2_%s(%s,xs,xp,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                           % ('SS', 'SS', all_3p_K2_ct[i][j], K2_3p_indices[j]))
                        list_str_M2_K2.append('if(ierr.eq.1)goto 999\n')
                        os.system('cat ' + NNLO_IR_limits_tmp_path + all_3p_K2_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
                    elif j == 1: # + SC_ijk
                        list_str_M2_K2.append('K%s=K%s+M2_%s(%s,iref,xs,xp,xsb,xpb,xsbb,xpbb,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                           % ('SC', 'SC', all_3p_K2_ct[i][j], K2_3p_indices[j]))
                        list_str_M2_K2.append('if(ierr.eq.1)goto 999\n')
                        os.system('cat ' + NNLO_IR_limits_tmp_path + all_3p_K2_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
                    elif j == 2: # - SS_ij SC_ijk
                        list_str_M2_K2.append('K%s=K%s-M2_%s(%s,iref,xs,xp,xsb,xpb,xsbb,xpbb,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                           % ('SC', 'SC', all_3p_K2_ct[i][j], K2_3p_indices[j]))
                        list_str_M2_K2.append('if(ierr.eq.1)goto 999\n')
                        os.system('cat ' + NNLO_IR_limits_tmp_path + all_3p_K2_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
                    elif j == 3: # + CC_ijk
                        list_str_M2_K2.append('K%s=K%s+M2_%s(%s,iref,xs,xp,xsb,xpb,xsbb,xpbb,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                           % ('CC', 'CC', all_3p_K2_ct[i][j], K2_3p_indices[j]))
                        list_str_M2_K2.append('if(ierr.eq.1)goto 999\n')
                        os.system('cat ' + NNLO_IR_limits_tmp_path + all_3p_K2_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
                    elif j == 4: # - SS_ij CC_ijk
                        list_str_M2_K2.append('K%s=K%s-M2_%s(%s,iref,xs,xp,xsb,xpb,xsbb,xpbb,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                           % ('CC', 'CC', all_3p_K2_ct[i][j], K2_3p_indices[j]))
                        list_str_M2_K2.append('if(ierr.eq.1)goto 999\n')
                        os.system('cat ' + NNLO_IR_limits_tmp_path + all_3p_K2_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
                    elif j == 5: # - SC_ijk CC_ijk
                        list_str_M2_K2.append('K%s=K%s-M2_%s(%s,iref,xs,xp,xsb,xpb,xsbb,xpbb,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                           % ('CC', 'CC', all_3p_K2_ct[i][j], K2_3p_indices[j]))
                        list_str_M2_K2.append('if(ierr.eq.1)goto 999\n')
                        os.system('cat ' + NNLO_IR_limits_tmp_path + all_3p_K2_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
                    elif j == 6: # + SS_ij SC_ijk CC_ijk
                        list_str_M2_K2.append('K%s=K%s+M2_%s(%s,iref,xs,xp,xsb,xpb,xsbb,xpbb,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                           % ('CC', 'CC', all_3p_K2_ct[i][j], K2_3p_indices[j]))
                        list_str_M2_K2.append('if(ierr.eq.1)goto 999\n')
                        os.system('cat ' + NNLO_IR_limits_tmp_path + all_3p_K2_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
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

            # L2_ijkj : 11 -> [Sik, SCijk, SCkij, SikSCijk, SikSCkij, Cijk, SikCijk, SCijkCijk, SikSCijkCijk, SCkijCijk, SikSCkijCijk]
            elif label == 'ijkj':
                tmp_str = """
c       KSS = SS_ik
c       KSC = (SC_ijk + SC_kij) (1 - SS_ik)
c       KCC = CC_ijk (1 - SS_ik) (1 - SC_ijk - SC_kij)"""
                list_str_defK2.append(tmp_str)
                for j in range(0, len(all_3p_K2_ct[i])):
                    if all_3p_K2_ct[i][j] ==  0:
                        continue
                    if j == 0: # SS_ik
                        list_str_M2_K2.append('K%s=K%s+M2_%s(%s,xs,xp,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                           % ('SS', 'SS', all_3p_K2_ct[i][j], K2_3p_indices[j]))
                        list_str_M2_K2.append('if(ierr.eq.1)goto 999\n')
                        os.system('cat ' + NNLO_IR_limits_tmp_path + all_3p_K2_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
                    elif j == 1: # + SC_ijk
                        list_str_M2_K2.append('K%s=K%s+M2_%s(%s,xs,xp,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                           % ('SC', 'SC', all_3p_K2_ct[i][j], K2_3p_indices[j]))
                        list_str_M2_K2.append('if(ierr.eq.1)goto 999\n')
                        os.system('cat ' + NNLO_IR_limits_tmp_path + all_3p_K2_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
                    elif j == 2: # + SC_kij
                        list_str_M2_K2.append('K%s=K%s+M2_%s(%s,xs,xp,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                           % ('SC', 'SC', all_3p_K2_ct[i][j], K2_3p_indices[j]))
                        list_str_M2_K2.append('if(ierr.eq.1)goto 999\n')
                        os.system('cat ' + NNLO_IR_limits_tmp_path + all_3p_K2_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
                    elif j == 3: # - SS_ik SC_ijk
                        list_str_M2_K2.append('K%s=K%s-M2_%s(%s,xs,xp,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                           % ('SC', 'SC', all_3p_K2_ct[i][j], K2_3p_indices[j]))
                        list_str_M2_K2.append('if(ierr.eq.1)goto 999\n')
                        os.system('cat ' + NNLO_IR_limits_tmp_path + all_3p_K2_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
                    elif j == 4: # - SS_ik SC_kij
                        list_str_M2_K2.append('K%s=K%s-M2_%s(%s,xs,xp,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                           % ('SC', 'SC', all_3p_K2_ct[i][j], K2_3p_indices[j]))
                        list_str_M2_K2.append('if(ierr.eq.1)goto 999\n')
                        os.system('cat ' + NNLO_IR_limits_tmp_path + all_3p_K2_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
                    elif j == 5: # + CC_ijk
                        list_str_M2_K2.append('K%s=K%s+M2_%s(%s,iref,xs,xp,xsb,xpb,xsbb,xpbb,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                           % ('CC', 'CC', all_3p_K2_ct[i][j], K2_3p_indices[j]))
                        list_str_M2_K2.append('if(ierr.eq.1)goto 999\n')
                        os.system('cat ' + NNLO_IR_limits_tmp_path + all_3p_K2_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
                    elif j == 6: # - SS_ik CC_ijk
                        list_str_M2_K2.append('K%s=K%s-M2_%s(%s,xs,xp,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                           % ('CC', 'CC', all_3p_K2_ct[i][j], K2_3p_indices[j]))
                        list_str_M2_K2.append('if(ierr.eq.1)goto 999\n')
                        os.system('cat ' + NNLO_IR_limits_tmp_path + all_3p_K2_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
                    elif j == 7: # - SC_ijk CC_ijk
                        list_str_M2_K2.append('K%s=K%s-M2_%s(%s,xs,xp,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                           % ('CC', 'CC', all_3p_K2_ct[i][j], K2_3p_indices[j]))
                        list_str_M2_K2.append('if(ierr.eq.1)goto 999\n')
                        os.system('cat ' + NNLO_IR_limits_tmp_path + all_3p_K2_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
                    elif j == 8: # + SS_ik SC_ijk CC_ijk
                        list_str_M2_K2.append('K%s=K%s+M2_%s(%s,xs,xp,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                           % ('CC', 'CC', all_3p_K2_ct[i][j], K2_3p_indices[j]))
                        list_str_M2_K2.append('if(ierr.eq.1)goto 999\n')
                        os.system('cat ' + NNLO_IR_limits_tmp_path + all_3p_K2_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
                    elif j == 9: # - SC_kij CC_ijk
                        list_str_M2_K2.append('K%s=K%s-M2_%s(%s,xs,xp,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                           % ('CC', 'CC', all_3p_K2_ct[i][j], K2_3p_indices[j]))
                        list_str_M2_K2.append('if(ierr.eq.1)goto 999\n')
                        os.system('cat ' + NNLO_IR_limits_tmp_path + all_3p_K2_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
                    elif j == 10: # + SS_ik SC_kij CC_ijk
                        list_str_M2_K2.append('K%s=K%s+M2_%s(%s,iref,xs,xp,xsb,xpb,xsbb,xpbb,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                           % ('CC', 'CC', all_3p_K2_ct[i][j], K2_3p_indices[j]))
                        list_str_M2_K2.append('if(ierr.eq.1)goto 999\n')
                        os.system('cat ' + NNLO_IR_limits_tmp_path + all_3p_K2_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
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

            # L12_ijjk : 13 -> [Si Sij, Si SCijk, Si SijSCijk, Si Cijk, Si SijCijk, Si SCijkCijk, Si SijSCijkCijk,
            #                   Cij Sij, Cij Cijk, Cij SijCijk, SiCij Sij, SiCij Cijk, SiCij SijCijk]

            ct_list = []
            if label == 'ijjk':
                tmp_str = """
c       KS_SS  = S_i SS_ij
c       KS_SC  = S_i SC_ijk (1 - SS_ij)
c       KS_CC  = S_i CC_ijk (1 - SS_ij) (1 - SC_ijk)
c       KHC_SS = C_ij (1-S_i) SS_ij
c       KHC_SC = 0
c       KHC_CC = C_ij (1-S_i) CC_ijk (1 - SS_ij)"""
                list_str_defK12.append(tmp_str)
                for j in range(0, len(all_3p_K12_ct[i])):
                    if all_3p_K12_ct[i][j] == 0:
                        continue
                    if j == 0:    # S_i SS_ij
                        list_str_M2_K12.append('K%s=K%s+M2_%s(%s,iref,xs,xp,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                       % ('S_SS', 'S_SS', all_3p_K12_ct[i][j], K12_3p_indices[j]))
                        list_str_M2_K12.append('if(ierr.eq.1)goto 999\n')
                        os.system('cat ' + NNLO_IR_limits_tmp_path + all_3p_K12_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
                    elif j == 1:  # + S_i SC_ijk
                        list_str_M2_K12.append('K%s=K%s+M2_%s(%s,iref,xs,xp,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                       % ('S_SC', 'S_SC', all_3p_K12_ct[i][j], K12_3p_indices[j]))
                        list_str_M2_K12.append('if(ierr.eq.1)goto 999\n')
                        os.system('cat ' + NNLO_IR_limits_tmp_path + all_3p_K12_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
                    elif j == 2:  # - S_i SS_ij SC_ijk
                        list_str_M2_K12.append('K%s=K%s-M2_%s(%s,iref,xs,xp,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                       % ('S_SC', 'S_SC', all_3p_K12_ct[i][j], K12_3p_indices[j]))
                        list_str_M2_K12.append('if(ierr.eq.1)goto 999\n')
                        os.system('cat ' + NNLO_IR_limits_tmp_path + all_3p_K12_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
                    elif j == 3:  # + S_i CC_ijk
                        list_str_M2_K12.append('K%s=K%s+M2_%s(%s,iref,xs,xp,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                       % ('S_CC', 'S_CC', all_3p_K12_ct[i][j], K12_3p_indices[j]))
                        list_str_M2_K12.append('if(ierr.eq.1)goto 999\n')
                        os.system('cat ' + NNLO_IR_limits_tmp_path + all_3p_K12_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
                    elif j == 4:  # - S_i SS_ij CC_ijk
                        list_str_M2_K12.append('K%s=K%s-M2_%s(%s,iref,xs,xp,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                       % ('S_CC', 'S_CC', all_3p_K12_ct[i][j], K12_3p_indices[j]))
                        list_str_M2_K12.append('if(ierr.eq.1)goto 999\n')
                        os.system('cat ' + NNLO_IR_limits_tmp_path + all_3p_K12_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
                    elif j == 5:  # - S_i SC_ijk CC_ijk
                        list_str_M2_K12.append('K%s=K%s-M2_%s(%s,iref,xs,xp,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                       % ('S_CC', 'S_CC', all_3p_K12_ct[i][j], K12_3p_indices[j]))
                        list_str_M2_K12.append('if(ierr.eq.1)goto 999\n')
                        os.system('cat ' + NNLO_IR_limits_tmp_path + all_3p_K12_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
                    elif j == 6:  # + S_i SS_ij SC_ijk CC_ijk
                        list_str_M2_K12.append('K%s=K%s+M2_%s(%s,iref,xs,xp,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                       % ('S_CC', 'S_CC', all_3p_K12_ct[i][j], K12_3p_indices[j]))
                        list_str_M2_K12.append('if(ierr.eq.1)goto 999\n')
                        os.system('cat ' + NNLO_IR_limits_tmp_path + all_3p_K12_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
                    elif j == 7:  # + C_ij SS_ij
                        list_str_M2_K12.append('K%s=K%s+M2_%s(%s,iref,xs,xp,xsb,xpb,xsbb,xpbb,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                       % ('HC_SS', 'HC_SS', all_3p_K12_ct[i][j], K12_3p_indices[j]))
                        list_str_M2_K12.append('if(ierr.eq.1)goto 999\n')
                        os.system('cat ' + NNLO_IR_limits_tmp_path + all_3p_K12_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
                    elif j == 8:  # + C_ij CC_ijk
                        list_str_M2_K12.append('K%s=K%s+M2_%s(%s,iref,xs,xp,xsb,xpb,xsbb,xpbb,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                       % ('HC_CC', 'HC_CC', all_3p_K12_ct[i][j], K12_3p_indices[j]))
                        list_str_M2_K12.append('if(ierr.eq.1)goto 999\n')
                        os.system('cat ' + NNLO_IR_limits_tmp_path + all_3p_K12_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
                    elif j == 9: # - C_ij SS_ij CC_ijk
                        list_str_M2_K12.append('K%s=K%s-M2_%s(%s,iref,xs,xp,xsb,xpb,xsbb,xpbb,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                       % ('HC_CC', 'HC_CC', all_3p_K12_ct[i][j], K12_3p_indices[j]))
                        list_str_M2_K12.append('if(ierr.eq.1)goto 999\n')
                        os.system('cat ' + NNLO_IR_limits_tmp_path + all_3p_K12_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
                    elif j == 10: # - S_i C_ij SS_ij
                        list_str_M2_K12.append('K%s=K%s-M2_%s(%s,iref,xs,xp,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                       % ('HC_SS', 'HC_SS', all_3p_K12_ct[i][j], K12_3p_indices[j]))
                        list_str_M2_K12.append('if(ierr.eq.1)goto 999\n')
                        os.system('cat ' + NNLO_IR_limits_tmp_path + all_3p_K12_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
                    elif j == 11: # - S_i C_ij CC_ijk
                        list_str_M2_K12.append('K%s=K%s-M2_%s(%s,iref,xs,xp,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                       % ('HC_CC', 'HC_CC', all_3p_K12_ct[i][j], K12_3p_indices[j]))
                        list_str_M2_K12.append('if(ierr.eq.1)goto 999\n')
                        os.system('cat ' + NNLO_IR_limits_tmp_path + all_3p_K12_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
                    elif j == 12: # + S_i C_ij SS_ij CC_ijk
                        list_str_M2_K12.append('K%s=K%s+M2_%s(%s,iref,xs,xp,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                       % ('HC_CC', 'HC_CC', all_3p_K12_ct[i][j], K12_3p_indices[j]))
                        list_str_M2_K12.append('if(ierr.eq.1)goto 999\n')
                        os.system('cat ' + NNLO_IR_limits_tmp_path + all_3p_K12_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')

                    if all_3p_K12_ct[i][j] not in ct_list:
                        ct_list.append(all_3p_K12_ct[i][j])
                        tmp_str = """
c       %s
        DOUBLE PRECISION M2_%s""" %(K12_labels[j],all_3p_K12_ct[i][j])
                        list_str_defK12.append(tmp_str)

            elif label == 'ijkj':
                # L12_ijkj : 13 -> [Si Sik, Si SCijk, Si SikSCijk, Si Cijk, Si SikCijk, Si SCijkCijk, Si SikSCijkCijk,
                #                   Cij SCkij, Cij Cijk, Cij SCkijCijk, SiCij SCkij, SiCij Cijk, SiCij SCkijCijk]
                tmp_str = """
c       KS_SS  = S_i SS_ik
c       KS_SC  = S_i SC_ijk (1 - SS_ik)
c       KS_CC  = S_i CC_ijk (1 - SS_ik) (1 - SC_ijk)
c       KHC_SS = 0
c       KHC_SC = C_ij (1-S_i) SC_kij
c       KHC_CC = C_ij (1-S_i) CC_ijk (1 - SC_kij)"""
                list_str_defK12.append(tmp_str)
                for j in range(0, len(all_3p_K12_ct[i])):
                    if all_3p_K12_ct[i][j] == 0:
                        continue
                    if j == 0:    # S_i SS_ik
                        list_str_M2_K12.append('K%s=K%s+M2_%s(isec,%s,iref,xs,xp,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                       % ('S_SS', 'S_SS', all_3p_K12_ct[i][j], K12_3p_indices[j]))
                        list_str_M2_K12.append('if(ierr.eq.1)goto 999\n')
                        os.system('cat ' + NNLO_IR_limits_tmp_path + all_3p_K12_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
                    elif j == 1:  # + S_i SC_ijk
                        list_str_M2_K12.append('K%s=K%s+M2_%s(isec,%s,iref,xs,xp,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                       % ('S_SC', 'S_SC', all_3p_K12_ct[i][j], K12_3p_indices[j]))
                        list_str_M2_K12.append('if(ierr.eq.1)goto 999\n')
                        os.system('cat ' + NNLO_IR_limits_tmp_path + all_3p_K12_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
                    elif j == 2:  # - S_i SS_ik SC_ijk
                        list_str_M2_K12.append('K%s=K%s-M2_%s(%s,iref,xs,xp,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                       % ('S_SC', 'S_SC', all_3p_K12_ct[i][j], K12_3p_indices[j]))
                        list_str_M2_K12.append('if(ierr.eq.1)goto 999\n')
                        os.system('cat ' + NNLO_IR_limits_tmp_path + all_3p_K12_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
                    elif j == 3:  # + S_i CC_ijk
                        list_str_M2_K12.append('K%s=K%s+M2_%s(%s,iref,xs,xp,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                       % ('S_CC', 'S_CC', all_3p_K12_ct[i][j], K12_3p_indices[j]))
                        list_str_M2_K12.append('if(ierr.eq.1)goto 999\n')
                        os.system('cat ' + NNLO_IR_limits_tmp_path + all_3p_K12_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
                    elif j == 4:  # - S_i SS_ik CC_ijk
                        list_str_M2_K12.append('K%s=K%s-M2_%s(%s,iref,xs,xp,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                       % ('S_CC', 'S_CC', all_3p_K12_ct[i][j], K12_3p_indices[j]))
                        list_str_M2_K12.append('if(ierr.eq.1)goto 999\n')
                        os.system('cat ' + NNLO_IR_limits_tmp_path + all_3p_K12_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
                    elif j == 5:  # - S_i SC_ijk CC_ijk
                        list_str_M2_K12.append('K%s=K%s-M2_%s(%s,iref,xs,xp,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                       % ('S_CC', 'S_CC', all_3p_K12_ct[i][j], K12_3p_indices[j]))
                        list_str_M2_K12.append('if(ierr.eq.1)goto 999\n')
                        os.system('cat ' + NNLO_IR_limits_tmp_path + all_3p_K12_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
                    elif j == 6:  # + S_i SS_ik SC_ijk CC_ijk
                        list_str_M2_K12.append('K%s=K%s+M2_%s(%s,iref,xs,xp,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                       % ('S_CC', 'S_CC', all_3p_K12_ct[i][j], K12_3p_indices[j]))
                        list_str_M2_K12.append('if(ierr.eq.1)goto 999\n')
                        os.system('cat ' + NNLO_IR_limits_tmp_path + all_3p_K12_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
                    elif j == 7:  # + C_ij SC_kij
                        list_str_M2_K12.append('K%s=K%s+M2_%s(%s,iref,xs,xp,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                       % ('HC_SC', 'HC_SC', all_3p_K12_ct[i][j], K12_3p_indices[j]))
                        list_str_M2_K12.append('if(ierr.eq.1)goto 999\n')
                        os.system('cat ' + NNLO_IR_limits_tmp_path + all_3p_K12_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
                    elif j == 8:  # + C_ij CC_ijk
                        list_str_M2_K12.append('K%s=K%s+M2_%s(%s,iref,xs,xp,xsb,xpb,xsbb,xpbb,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                       % ('HC_CC', 'HC_CC', all_3p_K12_ct[i][j], K12_3p_indices[j]))
                        list_str_M2_K12.append('if(ierr.eq.1)goto 999\n')
                        os.system('cat ' + NNLO_IR_limits_tmp_path + all_3p_K12_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
                    elif j == 9: # - C_ij SC_kij CC_ijk
                        list_str_M2_K12.append('K%s=K%s-M2_%s(%s,iref,xs,xp,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                       % ('HC_CC', 'HC_CC', all_3p_K12_ct[i][j], K12_3p_indices[j]))
                        list_str_M2_K12.append('if(ierr.eq.1)goto 999\n')
                        os.system('cat ' + NNLO_IR_limits_tmp_path + all_3p_K12_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
                    elif j == 10: # - S_i C_ij SC_kij
                        list_str_M2_K12.append('K%s=K%s-M2_%s(%s,iref,xs,xp,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                       % ('HC_SC', 'HC_SC', all_3p_K12_ct[i][j], K12_3p_indices[j]))
                        list_str_M2_K12.append('if(ierr.eq.1)goto 999\n')
                        os.system('cat ' + NNLO_IR_limits_tmp_path + all_3p_K12_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
                    elif j == 11: # - S_i C_ij CC_ijk
                        list_str_M2_K12.append('K%s=K%s-M2_%s(%s,iref,xs,xp,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                       % ('HC_CC', 'HC_CC', all_3p_K12_ct[i][j], K12_3p_indices[j]))
                        list_str_M2_K12.append('if(ierr.eq.1)goto 999\n')
                        os.system('cat ' + NNLO_IR_limits_tmp_path + all_3p_K12_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
                    elif j == 12: # + S_i C_ij SC_kij CC_ijk
                        list_str_M2_K12.append('K%s=K%s+M2_%s(%s,iref,xs,xp,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                       % ('HC_CC', 'HC_CC', all_3p_K12_ct[i][j], K12_3p_indices[j]))
                        list_str_M2_K12.append('if(ierr.eq.1)goto 999\n')
                        os.system('cat ' + NNLO_IR_limits_tmp_path + all_3p_K12_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')

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

            # write NNLO_K
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
            if label == 'ijjk':
                filename = pjoin(dirpath, 'NNLO_K_%d_%d_%d_%d.f' % (isec, jsec, jsec, ksec))
            elif label == 'ijkj':
                filename = pjoin(dirpath, 'NNLO_K_%d_%d_%d_%d.f' % (isec, jsec, ksec, jsec))
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
            # print(' ')

            # write NNLO_RRsub
            if sector_info['Born_str']:
                replace_dict_double_real['UBgraphs'] = overall_sector_info[i]['Born_str']
                replace_dict_limits['proc_prefix_Born'] = overall_sector_info[i]['Born_str']
            else:
                # Set dummy calls to bypass Born multichannelling
                if len(glob.glob("%s/ngraphs_dummy.inc" % dirpath)) == 0:
                    os.symlink(dirmadnklo + '/Template/Fortran_tmp/src_to_common/matrix_dummy.f', dirpath + '/matrix_dummy.f')
                    os.symlink(dirmadnklo + '/Template/Fortran_tmp/src_to_common/ngraphs_dummy.inc', dirpath + '/ngraphs_dummy.inc')
                    os.symlink(dirmadnklo + '/Template/Fortran_tmp/src_to_common/configs_dummy.f', dirpath + '/configs_dummy.f')
                    os.symlink(dirmadnklo + '/Template/Fortran_tmp/src_to_common/props_dummy.f', dirpath + '/props_dummy.f')
                    os.symlink(dirmadnklo + '/Template/Fortran_tmp/src_to_common/decayBW_dummy.f', dirpath + '/decayBW_dummy.f')
                    os.symlink(dirmadnklo + '/Template/Fortran_tmp/src_to_common/leshouche_dummy.f', dirpath + '/leshouche_dummy.f')
                    os.symlink(dirmadnklo + '/Template/Fortran_tmp/src_to_common/dummy_multich.f', dirpath + '/dummy_multich.f')

            filename = []
            if label == 'ijjk':
                filename = pjoin(dirpath, 'NNLO_RRsub_%d_%d_%d_%d.f' % (isec, jsec, jsec, ksec))
            elif label == 'ijkj':
                filename = pjoin(dirpath, 'NNLO_RRsub_%d_%d_%d_%d.f' % (isec, jsec, ksec, jsec))
            file = open(pjoin(dirmadnklo,"tmp_fortran/tmp_files/NNLO_RRsub_template.f")).read()
            file = file % replace_dict_double_real
            writer(filename).writelines(file)

            # write driver_RR
            UBgraphs = overall_sector_info[i]['Born_str']
            if label == 'ijjk':
                self.write_driver_npt_template(writer, dirpath, dirmadnklo, i , isec, jsec, jsec, ksec, UBgraphs)
            elif label == 'ijkj':
                self.write_driver_npt_template(writer, dirpath, dirmadnklo, i , isec, jsec, ksec, jsec, UBgraphs)

            # write test_RR
            if label == 'ijjk':
                self.write_testRR_template_file(writer, dirpath, dirmadnklo, defining_process,
                                   i, isec, jsec, jsec, ksec, all_3p_K1_ct, all_3p_K2_ct)
            elif label == 'ijkj':
                self.write_testRR_template_file(writer, dirpath, dirmadnklo, defining_process,
                                   i, isec, jsec, ksec, jsec, all_3p_K1_ct, all_3p_K2_ct)

            # write NNLO_IR_limits
            if label == 'ijjk':
                filename = pjoin(dirpath, 'NNLO_IR_limits_%d_%d_%d_%d.f' % (isec, jsec, jsec, ksec))
            elif label == 'ijkj':
                filename = pjoin(dirpath, 'NNLO_IR_limits_%d_%d_%d_%d.f' % (isec, jsec, ksec, jsec))
            file = open(NNLO_IR_limits_tmp_path + 'IR_tmp.f').read()
            file = file % replace_dict_limits
            writer(filename).writelines(file)
            # TODO: maybe safer to not remove this command
            os.system('rm ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')


######### Write NNLO_K_isec_jsec_ksec_lsec.f and NNLO_R_isec_jsec_ksec_lsec (4-particle sector)

        # Set replace_dict
        replace_dict_ct = {}
        replace_dict_limits = {}
        replace_dict_double_real ={}

        necessary_default_4p_ct_list = ['S_g', 'C_gg', 'C_gq', 'C_qqx', 'S_C_gg', 'S_C_gq', \
                                        'SS_gg', 'SS_qqx', \
                                        'SC_ggg', 'SC_ggq', 'SC_gqqx', \
                                        'SS_gg_SC_ggg', 'SS_gg_SC_ggq', \
                                        'CC_gggg', 'CC_gggq', 'CC_ggqqx', 'CC_qxqqxq', 'CC_gqgq', 'CC_gqqxq', \
                                        'CC_gqgg', 'CC_qqxgg', 'CC_qxqgq', \
                                        'SS_gg_CC_gggg', 'SS_gg_CC_gggq', 'SS_gg_CC_ggqqx', 'SS_qqx_CC_qxqqxq', \
                                        'SC_ggg_CC_gggg', 'SC_ggq_CC_gggq', 'SC_gqqx_CC_ggqqx', 'SC_ggq_CC_gqgq', 'SC_gqqx_CC_gqqxq', \
                                        'SC_ggq_CC_gqgg', 'SC_gqqx_CC_qqxgg', 'SC_gqqx_CC_qxqgq', \
                                        'S_SS_gg', 'S_SC_ggg', 'S_SC_ggq', 'S_SC_gqqx', \
                                        'S_SS_gg_SC_ggg', 'S_SS_gg_SC_ggq', \
                                        'C_SC_ggg', 'C_SC_ggq', 'C_SC_gqqx', \
                                        'C_CC_gggg', 'C_CC_gggq', 'C_CC_ggqqx', 'C_CC_qxqqxq', 'C_CC_gqgq', 'C_CC_gqqxq', \
                                        'C_CC_gqgg', 'C_CC_qqxgg', 'C_CC_qxqgq', \
                                        'C_SC_ggg_CC_gggg', 'C_SC_ggq_CC_gggq', 'C_SC_gqqx_CC_ggqqx', 'C_SC_ggq_CC_gqgq', 'C_SC_gqqx_CC_gqqxq', \
                                        'C_SC_ggq_CC_gqgg', 'C_SC_gqqx_CC_qqxgg', 'C_SC_gqqx_CC_qxqgq', \
                                        'SC_SC_ggg', 'SC_SC_ggq', 'SC_SC_gqqx', \
                                        'SC_CC_gggg', 'SC_CC_gggq', 'SC_CC_ggqqx', 'SC_CC_gqgq', 'SC_CC_gqqxq', 'SC_CC_gqgg', \
                                        'SC_SC_ggg_CC_gggg', 'SC_SC_ggq_CC_gggq', 'SC_SC_gqqx_CC_ggqqx', 'SC_SC_ggq_CC_gqgq', 'SC_SC_gqqx_CC_gqqxq', \
                                        'SC_SC_ggq_CC_gqgg']

        # K2_ijkl  : 9  ->  [Sik, SCikl, SCkij, SikSCikl, SikSCkij, Cijkl, SikCijkl, SCiklCijkl, SCkijCijkl]
        # K12_ijkl : 9  ->  [Si Sik, Si SCikl, Si SikSCikl, Cij SCkij, Cij Cijkl, Cij SCkijCijkl, SiCij SCkij, SiCij Cijkl, SiCij SCkijCijkl]
        K1_labels = ['S_i', 'C_ij', 'S_i C_ij']
        K2_labels = ['SS_ik', 'SC_ikl', 'SC_kij', 'SS_ik SC_ikl', 'SS_ik SC_kij', \
                     'CC_ijkl', 'SS_ik CC_ijkl', 'SC_ikl CC_ijkl', 'SC_kij CC_ijkl']
        K12_labels = ['S_i SS_ik', 'S_i SC_ikl', 'S_i SS_ik SC_ikl', \
                      'C_ij SC_kij', 'C_ij CC_ijkl', 'C_ij SC_kij CC_ijkl', \
                      'S_i C_ij SC_kij', 'S_i C_ij CC_ijkl', 'S_i C_ij SC_kij CC_ijkl']

        # Rule
        # CC_ijkl        -> isec,jsec,ksec,lsec
        # SS_ik CC_ijkl  -> isec,jsec,ksec,lsec
        # SC_ikl CC_ijkl -> isec,jsec,ksec,lsec
        # SC_kij CC_ijkl -> ksec,lsec,isec,jsec
        K1_4p_indices = ['isec', 'isec,jsec', 'isec,jsec']
        K2_4p_indices = ['isec,ksec', 'isec,ksec,lsec', 'ksec,isec,jsec', 'isec,ksec,lsec', 'ksec,isec,jsec', \
                         'isec,jsec,ksec,lsec', 'isec,jsec,ksec,lsec', 'isec,jsec,ksec,lsec', 'ksec,lsec,isec,jsec']
        K12_4p_indices = ['isec,ksec', 'isec,ksec,lsec', 'isec,ksec,lsec', \
                          'ksec,isec,jsec', 'isec,jsec,ksec,lsec', 'ksec,lsec,isec,jsec', \
                          'ksec,isec,jsec', 'isec,jsec,ksec,lsec', 'ksec,lsec,isec,jsec' ]

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
            replace_dict_ct['c3p'] = ksec
            replace_dict_ct['d3p'] = lsec
            replace_dict_double_real['isec'] = isec
            replace_dict_double_real['jsec'] = jsec
            replace_dict_double_real['ksec'] = ksec
            replace_dict_double_real['lsec'] = lsec
            replace_dict_double_real['iref'] = iref
            replace_dict_double_real['c3p'] = ksec
            replace_dict_double_real['d3p'] = lsec

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

            # Update sector_info dictionary
            sector_info = {
                'isec'          :   0,
                'jsec'          :   0,
                'ksec'          :   0,
                'lsec'          :   0,
                'iref'          :   0,
                'c3p'           :   0,
                'd3p'           :   0,
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
            sector_info['c3p'] = ksec
            sector_info['d3p'] = lsec

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

            # Initialise NLO_IR_limits.f for every sector [ij]
            list = []
            list.append('c Collection of relevant limits for sector [%d,%d,%d,%d] \n' %(isec,jsec,ksec,lsec))
            list.append('c K1  ' + str(all_4p_K1_ct[i])  + ' \n')
            list.append('c K2  ' + str(all_4p_K2_ct[i])  + ' \n')
            list.append('c K12 ' + str(all_4p_K12_ct[i]) + ' \n')
            string = "".join(list)
            NNLO_IR_limits_tmp_path = dirmadnklo + '/tmp_fortran/tmp_files/NNLO_limits/'
            with open(NNLO_IR_limits_tmp_path + "IR_tmp.f", "w") as f:
                f.writelines(string)

            # Loop on K1 cts
            # Recall that all_4p_K1_ct = [Si, Cij, SiCij]
            ct_list = []
            for j in range(0, len(all_4p_K1_ct[i])):
                if all_4p_K1_ct[i][j] ==  0:
                    continue
                if j == 0: # S_i
                    list_str_M2_K1.append('K%s=K%s+M2_%s(%s,xs,xp,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                       % (all_4p_K1_ct[i][j].split("_")[0], all_4p_K1_ct[i][j].split("_")[0], all_4p_K1_ct[i][j], K1_4p_indices[j]))
                    list_str_M2_K1.append('if(ierr.eq.1)goto 999\n')
                    os.system('cat ' + NNLO_IR_limits_tmp_path + all_4p_K1_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
                elif j == 1: # C_ij
                    list_str_M2_K1.append('K%s=K%s+M2_%s(%s,iref,xs,xp,xsb,xpb,wgt,xj,nitRR,1d0,wgt_chan,ierr)\n'
                                       % (all_4p_K1_ct[i][j].split("_")[0], all_4p_K1_ct[i][j].split("_")[0], all_4p_K1_ct[i][j], K1_4p_indices[j]))
                    list_str_M2_K1.append('if(ierr.eq.1)goto 999\n')
                    os.system('cat ' + NNLO_IR_limits_tmp_path + all_4p_K1_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
                elif j == 2: # S_C_ij
                    list_str_M2_K1.append('K%s=K%s+M2_%s(%s,xs,xp,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                       % ('SC', 'SC', all_4p_K1_ct[i][j], K1_4p_indices[j]))
                    list_str_M2_K1.append('if(ierr.eq.1)goto 999\n')
                    os.system('cat ' + NNLO_IR_limits_tmp_path + all_4p_K1_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')

                # Extract underlying real string
                self.get_uproc_str('Real', uB_all_4p_K1_ct[i][j], all_4p_K1_ct[i][j], dirpathR_head, replace_dict_limits,
                                       replace_dict_double_real, UReal_procs, path_UReal_procs, sector_info)

                if all_4p_K1_ct[i][j] not in ct_list:
                    ct_list.append(all_4p_K1_ct[i][j])
                    tmp_str = """
c       %s
        DOUBLE PRECISION M2_%s""" %(K1_labels[j],all_4p_K1_ct[i][j])
                    list_str_defK1.append(tmp_str)


            # Loop on K2 cts
            # K2_ijkl  : 9  ->  [Sik, SCikl, SCkij, SikSCikl, SikSCkij, Cijkl, SikCijkl, SCiklCijkl, SCkijCijkl]
            ct_list = []
            tmp_str = """
c       KSS = SS_ik
c       KSC = SC_ikl (1 - SS_ik) + SC_kij (1 - SS_ik)
c       KCC = CC_ijkl (1 + SS_ik - SC_ikl - SC_kij)"""
            list_str_defK2.append(tmp_str)
            for j in range(0, len(all_4p_K2_ct[i])):
                if all_4p_K2_ct[i][j] ==  0:
                    continue
                if j == 0: # SS_ik
                    list_str_M2_K2.append('K%s=K%s+M2_%s(%s,xs,xp,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                           % ('SS', 'SS', all_4p_K2_ct[i][j], K2_4p_indices[j]))
                    list_str_M2_K2.append('if(ierr.eq.1)goto 999\n')
                    os.system('cat ' + NNLO_IR_limits_tmp_path + all_4p_K2_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
                elif j == 1: # + SC_ikl
                    list_str_M2_K2.append('K%s=K%s+M2_%s(%s,xs,xp,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                           % ('SC', 'SC', all_4p_K2_ct[i][j], K2_4p_indices[j]))
                    list_str_M2_K2.append('if(ierr.eq.1)goto 999\n')
                    os.system('cat ' + NNLO_IR_limits_tmp_path + all_4p_K2_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
                elif j == 2: # + SC_kij
                    list_str_M2_K2.append('K%s=K%s+M2_%s(%s,xs,xp,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                           % ('SC', 'SC', all_4p_K2_ct[i][j], K2_4p_indices[j]))
                    list_str_M2_K2.append('if(ierr.eq.1)goto 999\n')
                    os.system('cat ' + NNLO_IR_limits_tmp_path + all_4p_K2_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
                elif j == 3: # - SS_ik SC_ikl
                    list_str_M2_K2.append('K%s=K%s-M2_%s(%s,xs,xp,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                           % ('SC', 'SC', all_4p_K2_ct[i][j], K2_4p_indices[j]))
                    list_str_M2_K2.append('if(ierr.eq.1)goto 999\n')
                    os.system('cat ' + NNLO_IR_limits_tmp_path + all_4p_K2_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
                elif j == 4: # - SS_ik SC_kij
                    list_str_M2_K2.append('K%s=K%s-M2_%s(%s,xs,xp,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                           % ('SC', 'SC', all_4p_K2_ct[i][j], K2_4p_indices[j]))
                    list_str_M2_K2.append('if(ierr.eq.1)goto 999\n')
                    os.system('cat ' + NNLO_IR_limits_tmp_path + all_4p_K2_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
                elif j == 5: # + CC_ijkl
                    list_str_M2_K2.append('K%s=K%s+M2_%s(%s,xs,xp,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                           % ('CC', 'CC', all_4p_K2_ct[i][j], K2_4p_indices[j]))
                    list_str_M2_K2.append('if(ierr.eq.1)goto 999\n')
                    os.system('cat ' + NNLO_IR_limits_tmp_path + all_4p_K2_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
                elif j == 6: # + SS_ik CC_ijkl
                    list_str_M2_K2.append('K%s=K%s+M2_%s(%s,xs,xp,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                           % ('CC', 'CC', all_4p_K2_ct[i][j], K2_4p_indices[j]))
                    list_str_M2_K2.append('if(ierr.eq.1)goto 999\n')
                    os.system('cat ' + NNLO_IR_limits_tmp_path + all_4p_K2_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
                elif j == 7: # - SC_ikl CC_ijkl
                    list_str_M2_K2.append('K%s=K%s-M2_%s(%s,xs,xp,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                           % ('CC', 'CC', all_4p_K2_ct[i][j], K2_4p_indices[j]))
                    list_str_M2_K2.append('if(ierr.eq.1)goto 999\n')
                    os.system('cat ' + NNLO_IR_limits_tmp_path + all_4p_K2_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
                elif j == 8: # - SC_kij CC_ijkl
                    list_str_M2_K2.append('K%s=K%s-M2_%s(%s,xs,xp,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                           % ('CC', 'CC', all_4p_K2_ct[i][j], K2_4p_indices[j]))
                    list_str_M2_K2.append('if(ierr.eq.1)goto 999\n')
                    os.system('cat ' + NNLO_IR_limits_tmp_path + all_4p_K2_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
                # Extract underlying Born string
                self.get_uproc_str('Born', uB_all_4p_K2_ct[i][j], all_4p_K2_ct[i][j], dirpathB_head, replace_dict_limits,
                                       replace_dict_double_real, UBorn_procs, path_UBorn_procs, sector_info)

                if all_4p_K2_ct[i][j] not in ct_list:
                    ct_list.append(all_4p_K2_ct[i][j])
                    tmp_str = """
c       %s
        DOUBLE PRECISION M2_%s""" %(K2_labels[j],all_4p_K2_ct[i][j])
                    list_str_defK2.append(tmp_str)

            # Loop on K12 cts
            # K12_ijkl : 9  ->  [Si Sik, Si SCikl, Si SikSCikl, Cij SCkij, Cij Cijkl, Cij SCkijCijkl, SiCij SCkij, SiCij Cijkl, SiCij SCkijCijkl]
            ct_list = []
            tmp_str = """
c       KS_SS  = S_i SS_ik
c       KS_SC  = S_i SC_ikl (1 - SS_ik)
c       KS_CC  = 0
c       KHC_SS = 0
c       KHC_SC = C_ij (1-S_i) SC_kij
c       KHC_CC = C_ij (1-S_i) CC_ijkl (1 - SC_kij)"""
            list_str_defK12.append(tmp_str)
            for j in range(0, len(all_4p_K12_ct[i])):
                if all_4p_K12_ct[i][j] == 0:
                    continue
                if j == 0:    # S_i SS_ik
                    list_str_M2_K12.append('K%s=K%s+M2_%s(%s,iref,xs,xp,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                       % ('S_SS', 'S_SS', all_4p_K12_ct[i][j], K12_4p_indices[j]))
                    list_str_M2_K12.append('if(ierr.eq.1)goto 999\n')
                    os.system('cat ' + NNLO_IR_limits_tmp_path + all_4p_K12_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
                elif j == 1:  # + S_i SC_ikl
                    list_str_M2_K12.append('K%s=K%s+M2_%s(%s,iref,xs,xp,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                       % ('S_SC', 'S_SC', all_4p_K12_ct[i][j], K12_4p_indices[j]))
                    list_str_M2_K12.append('if(ierr.eq.1)goto 999\n')
                    os.system('cat ' + NNLO_IR_limits_tmp_path + all_4p_K12_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
                elif j == 2:  # - S_i SS_ik SC_ikl
                    list_str_M2_K12.append('K%s=K%s-M2_%s(%s,iref,xs,xp,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                       % ('S_SC', 'S_SC', all_4p_K12_ct[i][j], K12_4p_indices[j]))
                    list_str_M2_K12.append('if(ierr.eq.1)goto 999\n')
                    os.system('cat ' + NNLO_IR_limits_tmp_path + all_4p_K12_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
                elif j == 3:  # + C_ij SC_kij
                    list_str_M2_K12.append('K%s=K%s+M2_%s(%s,iref,xs,xp,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                       % ('HC_SC', 'HC_SC', all_4p_K12_ct[i][j], K12_4p_indices[j]))
                    list_str_M2_K12.append('if(ierr.eq.1)goto 999\n')
                    os.system('cat ' + NNLO_IR_limits_tmp_path + all_4p_K12_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
                elif j == 4:  # + C_ij CC_ijkl
                    list_str_M2_K12.append('K%s=K%s+M2_%s(%s,iref,xs,xp,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                       % ('HC_CC', 'HC_CC', all_4p_K12_ct[i][j], K12_4p_indices[j]))
                    list_str_M2_K12.append('if(ierr.eq.1)goto 999\n')
                    os.system('cat ' + NNLO_IR_limits_tmp_path + all_4p_K12_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
                elif j == 5:  # - C_ij SC_kij CC_ijkl
                    list_str_M2_K12.append('K%s=K%s-M2_%s(%s,iref,xs,xp,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                       % ('HC_CC', 'HC_CC', all_4p_K12_ct[i][j], K12_4p_indices[j]))
                    list_str_M2_K12.append('if(ierr.eq.1)goto 999\n')
                    os.system('cat ' + NNLO_IR_limits_tmp_path + all_4p_K12_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
                elif j == 6:  # - S_i C_ij SC_kij
                    list_str_M2_K12.append('K%s=K%s-M2_%s(%s,iref,xs,xp,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                       % ('HC_SC', 'HC_SC', all_4p_K12_ct[i][j], K12_4p_indices[j]))
                    list_str_M2_K12.append('if(ierr.eq.1)goto 999\n')
                    os.system('cat ' + NNLO_IR_limits_tmp_path + all_4p_K12_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
                elif j == 7:  # - S_i C_ij CC_ijkl
                    list_str_M2_K12.append('K%s=K%s-M2_%s(%s,iref,xs,xp,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                       % ('HC_CC', 'HC_CC', all_4p_K12_ct[i][j], K12_4p_indices[j]))
                    list_str_M2_K12.append('if(ierr.eq.1)goto 999\n')
                    os.system('cat ' + NNLO_IR_limits_tmp_path + all_4p_K12_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')
                elif j == 8:  # + S_i C_ij SC_kij CC_ijkl
                    list_str_M2_K12.append('K%s=K%s+M2_%s(%s,iref,xs,xp,wgt,xj,xjB,nitRR,1d0,wgt_chan,ierr)\n'
                                       % ('HC_CC', 'HC_CC', all_4p_K12_ct[i][j], K12_4p_indices[j]))
                    list_str_M2_K12.append('if(ierr.eq.1)goto 999\n')
                    os.system('cat ' + NNLO_IR_limits_tmp_path + all_4p_K12_ct[i][j] + '.f >> ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')

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

            replace_dict_double_real['mapping_str'] = mapping_str

            # write NNLO_K
            filename = pjoin(dirpath, 'NNLO_K_%d_%d_%d_%d.f' % (isec, jsec, ksec, lsec))
            file = open(pjoin(dirmadnklo,"tmp_fortran/tmp_files/NNLO_K_template.f")).read()
            file = file % replace_dict_ct
            writer(filename).writelines(file)

#             # # check on sector_info
#             # print('Born_str : ' + str(overall_sector_info[i+len(all_3p_sector_list)]['Born_str']))
#             # print('alt_Born_str : ' + str(overall_sector_info[i+len(all_3p_sector_list)]['alt_Born_str']))
#             # print('Born_PDGs : ' + str(overall_sector_info[i+len(all_3p_sector_list)]['Born_PDGs']))
#             # print('path_to_Born : ' + str(overall_sector_info[i+len(all_3p_sector_list)]['path_to_Born']))
#             # print('alt_Born_path : ' + str(overall_sector_info[i+len(all_3p_sector_list)]['alt_Born_path']))
#             # print('Real_str : ' + str(overall_sector_info[i+len(all_3p_sector_list)]['Real_str']))
#             # print('Real_PDGs : ' + str(overall_sector_info[i+len(all_3p_sector_list)]['Real_PDGs']))
#             # print('path_to_Real : ' + str(overall_sector_info[i+len(all_3p_sector_list)]['path_to_Real']))
#             # print('alt_Real_str : ' + str(overall_sector_info[i+len(all_3p_sector_list)]['alt_Real_str']))
#             # print('alt_Real_path : ' + str(overall_sector_info[i+len(all_3p_sector_list)]['alt_Real_path']))
#             # print('  ')

            # write NNLO_RRsub
            if sector_info['Born_str']:
                replace_dict_double_real['UBgraphs'] = overall_sector_info[i]['Born_str']
                replace_dict_limits['proc_prefix_Born'] = overall_sector_info[i]['Born_str']
            else:
                # Set dummy calls to bypass Born multichannelling
                if len(glob.glob("%s/ngraphs_dummy.inc" % dirpath)) == 0:
                    os.symlink(dirmadnklo + '/Template/Fortran_tmp/src_to_common/matrix_dummy.f', dirpath + '/matrix_dummy.f')
                    os.symlink(dirmadnklo + '/Template/Fortran_tmp/src_to_common/ngraphs_dummy.inc', dirpath + '/ngraphs_dummy.inc')
                    os.symlink(dirmadnklo + '/Template/Fortran_tmp/src_to_common/configs_dummy.f', dirpath + '/configs_dummy.f')
                    os.symlink(dirmadnklo + '/Template/Fortran_tmp/src_to_common/props_dummy.f', dirpath + '/props_dummy.f')
                    os.symlink(dirmadnklo + '/Template/Fortran_tmp/src_to_common/decayBW_dummy.f', dirpath + '/decayBW_dummy.f')
                    os.symlink(dirmadnklo + '/Template/Fortran_tmp/src_to_common/leshouche_dummy.f', dirpath + '/leshouche_dummy.f')
                    os.symlink(dirmadnklo + '/Template/Fortran_tmp/src_to_common/dummy_multich.f', dirpath + '/dummy_multich.f')

            filename = []
            filename = pjoin(dirpath, 'NNLO_RRsub_%d_%d_%d_%d.f' % (isec, jsec, ksec, lsec))
            file = open(pjoin(dirmadnklo,"tmp_fortran/tmp_files/NNLO_RRsub_template.f")).read()
            file = file % replace_dict_double_real
            writer(filename).writelines(file)

            # write driver_RR
            UBgraphs = overall_sector_info[i+len(all_4p_sector_list)]['Born_str']
            self.write_driver_npt_template(writer, dirpath, dirmadnklo, i , isec, jsec, ksec, lsec, UBgraphs)

            # write test_RR
            self.write_testRR_template_file(writer, dirpath, dirmadnklo, defining_process,
                                   i, isec, jsec, ksec, lsec, all_4p_K1_ct, all_4p_K2_ct)

            # write NNLO_IR_limits
            filename = pjoin(dirpath, 'NNLO_IR_limits_%d_%d_%d_%d.f' % (isec, jsec, ksec, lsec))
            #file = open(pjoin(dirmadnklo,"tmp_fortran/tmp_files/NLO_IR_limits_tmp.f")).read()
            file = open(NNLO_IR_limits_tmp_path + 'IR_tmp.f').read()
            file = file % replace_dict_limits
            writer(filename).writelines(file)
            # TODO: maybe safer to reinstate this command
            os.system('rm ' + NNLO_IR_limits_tmp_path + 'IR_tmp.f')


# #---------- Functions outside loop on sectors ----------#
# ######### Write get_Born_PDGs.f & get_Real_PDGs.f

        self.write_get_Born_PDGs_file(writer, dirpath, overall_sector_info)
        self.write_get_Real_PDGs_file(writer, dirpath, overall_sector_info)
        self.write_get_UnderLying_PDGs_file(writer, dirpath, overall_sector_info)

######### Write makefile_npt_template

        self.write_makefile_RR_file(writer, dirpath, dirmadnklo, defining_process, overall_sector_info)

######### Write ajob_isec_jsec_ksec_lsec
        # Here stuff for multichanneling??


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
        # Write a DATA list including all sectors. Each sector is defined by 4 leg indices.

        # Define 3p sectors as 4p sectors with lsec = 0
        # all_3p_sector_list_with_0 = []
        # for i in range(0,len(all_3p_sector_list)):
        #     new_list = list(all_3p_sector_list[i])
        #     #new_list.append(0)
        #     new_list = tuple(new_list)
        #     all_3p_sector_list_with_0.append(new_list)
        # # Check
        # if len(all_3p_sector_list_with_0) != len(all_3p_sector_list):
        #     print('WARNING: Wrong number of 3p sectors!')
        #     return

        #all_sector_list = all_3p_sector_list_with_0 + all_4p_sector_list
        all_sector_list = all_3p_sector_list + all_4p_sector_list
        tmp_str = str(all_sector_list).replace('[','').replace(']','').replace(' ','').replace('(','').replace(')','')

        file = """ \
          integer, parameter :: lensectors = %d
          integer tmp(4,lensectors)
          integer all_sector_list(4,lensectors)""" % (len(all_sector_list))
        for k in range(0,len(all_sector_list)):
            file += """
          data tmp(1,%d), tmp(2,%d), tmp(3,%d), tmp(4,%d) /%s/ """ % (k+1,k+1,k+1,k+1,tmp_str[k*8:k*8+7])
          #  file += """
          #data all_sector_list(1,%d), all_sector_list(2,%d), all_sector_list(3,%d), all_sector_list(4,%d) /%s/ """ % (k+1,k+1,k+1,k+1,tmp_str[k*8:k*8+8])

        file += """
          all_sector_list = tmp"""

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
          integer Born_leg_PDGs(nexternal_UUB), Real_leg_PDGs(nexternal_UB)
          \n"""

        for i in range(0,len(overall_sector_info)):

            replace_dict_tmp = {}
            replace_dict_tmp['isec'] = overall_sector_info[i]['isec']
            replace_dict_tmp['jsec'] = overall_sector_info[i]['jsec']
            replace_dict_tmp['ksec'] = overall_sector_info[i]['c3p']
            replace_dict_tmp['lsec'] = overall_sector_info[i]['d3p']
            replace_dict_tmp['tmp_Real_PDGs'] = overall_sector_info[i]['Real_PDGs']
            replace_dict_tmp['tmp_Born_PDGs'] = overall_sector_info[i]['Born_PDGs']
            if(overall_sector_info[i]['Born_PDGs'] == []):
                replace_dict_tmp['tmp_Born_PDGs'] = '0'

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
          if(npart .eq. nexternal_UB) then
            Underlying_leg_PDGs = Real_leg_PDGs
          elseif(npart .eq. nexternal_UUB) then
            Underlying_leg_PDGs = Born_leg_PDGs
          else
            write(*,*) 'Get_Underlying_PDGs: error'
            write(*,*) 'npart is neither equal to nexternal_UB nor to nexternal_UUB:'
            write(*,*) 'npart, nexternal_UB, nexternal_UUB = ',  npart, nexternal_UB, nexternal_UUB
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

            # copy Born stuff
            if (not overall_sector_info[i]['Born_str'] and not overall_sector_info[i]['path_to_Born']):
                # No underlying Born
                continue

            elif (not overall_sector_info[i]['path_to_Born'] and overall_sector_info[i]['Born_str']):
                string = overall_sector_info[i]['Born_str']
                string2 = overall_sector_info[i]['alt_Born_str']
                # Set up link to matrix elements and their spin_correlation files related to the the flavour-dependent Born string
                if not glob.glob("%s/matrix_%s.f" % (dirpath, string2)):
                    os.symlink( "%s/matrix_%s.f" % (overall_sector_info[i]['alt_Born_path'], string2), "%s/matrix_%s.f" % (dirpath, string2) )
                    os.symlink( overall_sector_info[i]['alt_Born_path'] + '/%s_spin_correlations.inc' % string2, dirpath + '/%s_spin_correlations.inc' % string2 )

                if not glob.glob(dirpath + '/ngraphs_%s.inc' % string):
                    os.symlink(dirpath + '/../../../Common_Files/ngraphs_%s.inc' % string, dirpath + '/ngraphs_%s.inc' % string)
                    os.symlink(dirpath + '/../../../Common_Files/configs_%s.f' % string, dirpath + '/configs_%s.f' % string)
                    os.symlink(dirpath + '/../../../Common_Files/props_%s.f' % string, dirpath + '/props_%s.f' % string)
                    os.symlink(dirpath + '/../../../Common_Files/decayBW_%s.f' % string, dirpath + '/decayBW_%s.f' % string)
                    os.symlink(dirpath + '/../../../Common_Files/leshouche_%s.f' % string, dirpath + '/leshouche_%s.f' % string)

            elif (overall_sector_info[i]['path_to_Born'] and overall_sector_info[i]['Born_str']):
                string = overall_sector_info[i]['Born_str']
                # Set up link to matrix elements and their spin_correlation files related to the the flavour-dependent Born string
                if not glob.glob("%s/matrix_%s.f" % (dirpath, overall_sector_info[i]['Born_str'])):
                    os.symlink( "%s/matrix_%s.f" % (overall_sector_info[i]['path_to_Born'], string), "%s/matrix_%s.f" % (dirpath, string) )
                    os.symlink( overall_sector_info[i]['path_to_Born'] + '/%s_spin_correlations.inc' % string, dirpath + '/%s_spin_correlations.inc' % string )

                if not glob.glob(dirpath + '/ngraphs_%s.inc' % string):
                    os.symlink(dirpath + '/../../../Common_Files/ngraphs_%s.inc' % string, dirpath + '/ngraphs_%s.inc' % string)
                    os.symlink(dirpath + '/../../../Common_Files/configs_%s.f' % string, dirpath + '/configs_%s.f' % string)
                    os.symlink(dirpath + '/../../../Common_Files/props_%s.f' % string, dirpath + '/props_%s.f' % string)
                    os.symlink(dirpath + '/../../../Common_Files/decayBW_%s.f' % string, dirpath + '/decayBW_%s.f' % string)
                    os.symlink(dirpath + '/../../../Common_Files/leshouche_%s.f' % string, dirpath + '/leshouche_%s.f' % string)

            # copy Real stuff
            if (not overall_sector_info[i]['Real_str'] and not overall_sector_info[i]['path_to_Real']):
                # No underlying Real
                continue

            elif (not overall_sector_info[i]['path_to_Real'] and overall_sector_info[i]['Real_str']):
                string = overall_sector_info[i]['Real_str']
                string2 = overall_sector_info[i]['alt_Real_str']
                # Set up link to matrix elements and their spin_correlation files related to the the flavour-dependent Born string
                if not glob.glob("%s/matrix_%s.f" % (dirpath, string2)):
                    link_real = "%s/all_sector_list_real.inc" % dirpath
                    if os.path.lexists(link_real):
                        os.remove(link_real)
                    os.symlink( "%s/all_sector_list.inc" % overall_sector_info[i]['alt_Real_path'], link_real )
                    os.symlink( "%s/matrix_%s.f" % (overall_sector_info[i]['alt_Real_path'], string2), "%s/matrix_%s.f" % (dirpath, string2) )
                    os.symlink( overall_sector_info[i]['alt_Real_path'] + '/%s_spin_correlations.inc' % string2, dirpath + '/%s_spin_correlations.inc' % string2 )

                # if not glob.glob(dirpath + '/ngraphs_%s.inc' % string):
                #     os.symlink(dirpath + '/../../../Common_Files/ngraphs_%s.inc' % string, dirpath + '/ngraphs_%s.inc' % string)
                #     os.symlink(dirpath + '/../../../Common_Files/configs_%s.f' % string, dirpath + '/configs_%s.f' % string)
                #     os.symlink(dirpath + '/../../../Common_Files/props_%s.f' % string, dirpath + '/props_%s.f' % string)
                #     os.symlink(dirpath + '/../../../Common_Files/decayBW_%s.f' % string, dirpath + '/decayBW_%s.f' % string)
                #     os.symlink(dirpath + '/../../../Common_Files/leshouche_%s.f' % string, dirpath + '/leshouche_%s.f' % string)

            elif (overall_sector_info[i]['path_to_Real'] and overall_sector_info[i]['Real_str']):
                string = overall_sector_info[i]['Real_str']
                # Set up link to matrix elements and their spin_correlation files related to the the flavour-dependent Born string
                if not glob.glob("%s/matrix_%s.f" % (dirpath, overall_sector_info[i]['Real_str'])):
                    link_real = "%s/all_sector_list_real.inc" % dirpath
                    if os.path.lexists(link_real):
                        os.remove(link_real)
                    os.symlink( "%s/all_sector_list.inc" % overall_sector_info[i]['path_to_Real'], link_real )
                    os.symlink( "%s/matrix_%s.f" % (overall_sector_info[i]['path_to_Real'], string), "%s/matrix_%s.f" % (dirpath, string) )
                    os.symlink( overall_sector_info[i]['path_to_Real'] + '/%s_spin_correlations.inc' % string, dirpath + '/%s_spin_correlations.inc' % string )

                # if not glob.glob(dirpath + '/ngraphs_%s.inc' % string):
                #     os.symlink(dirpath + '/../../../Common_Files/ngraphs_%s.inc' % string, dirpath + '/ngraphs_%s.inc' % string)
                #     os.symlink(dirpath + '/../../../Common_Files/configs_%s.f' % string, dirpath + '/configs_%s.f' % string)
                #     os.symlink(dirpath + '/../../../Common_Files/props_%s.f' % string, dirpath + '/props_%s.f' % string)
                #     os.symlink(dirpath + '/../../../Common_Files/decayBW_%s.f' % string, dirpath + '/decayBW_%s.f' % string)
                #     os.symlink(dirpath + '/../../../Common_Files/leshouche_%s.f' % string, dirpath + '/leshouche_%s.f' % string)


    #===========================================================================
    # write file for testing limits, 'testRR.f'
    #===========================================================================

    def write_testRR_template_file(self, writer, dirpath, dirmadnklo, defining_process,
                                    i, isec, jsec, c3p, d3p, K1_ct, K2_ct):
        replace_dict = {}
        replace_dict['isec'] = isec
        replace_dict['jsec'] = jsec
        replace_dict['c3p'] = c3p
        replace_dict['d3p'] = d3p
        mass_list = []
        mass_list = 'ZERO' #FIRST ATTEMPT TO SKIP THE IF BELOW

        limit_str = ''

        # Mapping info
        if (jsec == c3p and jsec != d3p):     # ijjk: Sij, Cijk, SCijk
            limit_str += """
c
c     mapping ((i,j,r),(j,k,r))
"""
        elif (jsec != c3p and jsec == d3p):   # ijkj: Sik, Cijk, SCijk, SCkij
            limit_str += """
c
c     mapping ((i,j,r),(j,k,r))
"""
        elif (jsec != c3p and jsec != d3p):   # ijkl: Sik, Cijkl, SCikl, SCkij
            limit_str += """
c
c     mapping ((i,j,r),(k,l,r))
"""

        ### K1 limits to test: Si, Cij - common to ijjk,ijkj,ijkl ###
        # From (3.5) in 2212.11190
        if K1_ct[i][0] != 0 : #Si
            limit_str += """
c
c     limit Si
      e = [1d0,1d0,0d0,0d0,0d0]
      l = [0d0,0d0,0d0,0d0,0d0]
      call do_limit_RR_%d_%d_%d_%d(iunit,'Si      ',x0,e,l)
"""%(isec,jsec,c3p,d3p)

        # Loop over sectors with final state particles only
        if isec > 2 and jsec > 2:
            if K1_ct[i][1] != 0 : # Cij
                limit_str += """
c
c     limit Cij
      e=[0d0,1d0,0d0,0d0,0d0]
      l=[0d0,0d0,0d0,0d0,0d0]
      call do_limit_RR_%d_%d_%d_%d(iunit,'Cij     ',x0,e,l)
"""%(isec,jsec,c3p,d3p)

        elif isec > 2 and jsec <= 2:
            limit_str += """Collinear limits still to be specified in sectorsRR.py """
            raise MadEvent7Error('Collinear limits still to be specified in sectorsRR.py. ')

        ### K2 limits to test ###
        # K2_ijjk  : 7 ->  [Sij, SCijk, SijSCijk, Cijk, SijCijk, SCijkCijk, SijSCijkCijk]
        # K2_ijkj  : 11 -> [Sik, SCijk, SCkij, SikSCijk, SikSCkij, Cijk, SikCijk, SCijkCijk, SikSCijkCijk, SCkijCijk, SikSCkijCijk]
        # K2_ijkl  : 9  -> [Sik, SCikl, SCkij, SikSCikl, SikSCkij, Cijkl, SikCijkl, SCiklCijkl, SCkijCijkl]

        if (jsec == c3p and jsec != d3p):     # ijjk: Sij, SCijk, Cijk
            if K2_ct[i][0] != 0: # Sij
                limit_str += """
c
c     limit Sij
      e = [0d0,1d0,0d0,1d0,1d0]
      l = [0d0,0d0,0d0,0d0,0d0]
      call do_limit_RR_%d_%d_%d_%d(iunit,'Sij     ',x0,e,l)
"""%(isec,jsec,c3p,d3p)
            if K2_ct[i][1] != 0: # SCijk
                limit_str += """
c
c     limit SCijk
      e=[1d0,1d0,0d0,0d0,1d0]
      l=[0d0,0d0,0d0,0d0,0d0]
      call do_limit_RR_%d_%d_%d_%d(iunit,'SCijk   ',x0,e,l)
"""%(isec,jsec,c3p,d3p)
            if K2_ct[i][3] != 0: # Cijk
                limit_str += """
c
c       limit Cijk
        e=[0d0,1d0,0d0,0d0,1d0]
        l=[0d0,0d0,0d0,0d0,0d0]
      call do_limit_RR_%d_%d_%d_%d(iunit,'Cijk    ',x0,e,l)
"""%(isec,jsec,c3p,d3p)
        elif (jsec != c3p and jsec == d3p):   # ijkj: Sik, SCijk, SCkij, Cijk
            if K2_ct[i][0] != 0: # Sik
                limit_str += """
c
c     limit Sik
      e = [1d0,1d0,0d0,1d0,1d0]
      l = [0d0,0d0,0d0,1d0,0d0]
      call do_limit_RR_%d_%d_%d_%d(iunit,'Sik     ',x0,e,l)
"""%(isec,jsec,c3p,d3p)
            if K2_ct[i][1] != 0: # SCijk
                limit_str += """
c
c     limit SCijk
      e=[1d0,1d0,0d0,0d0,1d0]
      l=[0d0,0d0,0d0,0d0,0d0]
      call do_limit_RR_%d_%d_%d_%d(iunit,'SCijk   ',x0,e,l)
"""%(isec,jsec,c3p,d3p)
            if K2_ct[i][2] != 0: # SCkij
                limit_str += """
c
c     single soft double collinear limit SCkij
      e=[0d0,1d0,0d0,1d0,1d0]
      l=[0d0,0d0,0d0,1d0,0d0]
      call do_limit_RR_%d_%d_%d_%d(iunit,'SCkij   ',x0,e,l)
"""%(isec,jsec,c3p,d3p)
            if K2_ct[i][5] != 0: # Cijk
                limit_str += """
c
c       limit Cijk
        e=[0d0,1d0,0d0,0d0,1d0]
        l=[0d0,0d0,0d0,0d0,0d0]
      call do_limit_RR_%d_%d_%d_%d(iunit,'Cijk    ',x0,e,l)
"""%(isec,jsec,c3p,d3p)
        elif (jsec != c3p and jsec != d3p):   # ijkl: Sik, Cijkl, SCikl, SCkij
            limit_str += """
c       TODO: Testing limits for ijkl sector still to be specified in sectorsRR.py """



        # Test for spurious singularities
        # ijjk: 3 -> [Cir, Cjr, Cijr]
        # ijkj: 3 -> [Cir, Ckr, Cikr]
        # ijkl: 2 -> [Cir, Ckr]
        if (jsec == c3p and jsec != d3p):     # ijjk
            # Cir
            limit_str += """
c
c     spurious limit Cir
      e = [1d0,0d0,0d0,0d0,0d0]
      l = [0d0,0d0,0d0,0d0,0d0]
      call do_limit_RR_%d_%d_%d_%d(iunit,'Cir     ',x0,e,l)
"""%(isec,jsec,c3p,d3p)
            # Cjr
            limit_str += """
c
c     spurious limit Cjr
      e = [1d0,0d0,0d0,0d0,0d0]
      l = [1d0,0d0,0d0,0d0,0d0]
      call do_limit_RR_%d_%d_%d_%d(iunit,'Cjr     ',x0,e,l)
"""%(isec,jsec,c3p,d3p)
            # Cijr
            limit_str += """
c
c     spurious limit Cijr
      e = [0d0,0d0,0d0,1d0,0d0]
      l = [0d0,0d0,0d0,0d0,0d0]
      call do_limit_RR_%d_%d_%d_%d(iunit,'Cijr    ',x0,e,l)
"""%(isec,jsec,c3p,d3p)
        elif (jsec != c3p and jsec == d3p):   # ijkj
            # Cir
            limit_str += """
c
c     spurious limit Cir
      e = [1d0,0d0,0d0,0d0,0d0]
      l = [0d0,0d0,0d0,0d0,0d0]
      call do_limit_RR_%d_%d_%d_%d(iunit,'Cir     ',x0,e,l)
"""%(isec,jsec,c3p,d3p)
            # Ckr
            limit_str += """
c
c     spurious limit Ckr
      e = [0d0,0d0,0d0,1d0,0d0]
      l = [0d0,0d0,0d0,1d0,0d0]
      call do_limit_RR_%d_%d_%d_%d(iunit,'Ckr     ',x0,e,l)
"""%(isec,jsec,c3p,d3p)
            # Cikr
            limit_str += """
c
c     spurious limit Cikr
      e = [1d0,0d0,0d0,1d0,0d0]
      l = [0d0,0d0,0d0,1d0,0d0]
      call do_limit_RR_%d_%d_%d_%d(iunit,'Cikr    ',x0,e,l)
"""%(isec,jsec,c3p,d3p)
        elif (jsec != c3p and jsec != d3p):   # ijkl
            limit_str += """Testing limits for ijkl sector still to be specified in sectorsRR.py """


        replace_dict['limit_str'] = limit_str
        replace_dict['NNLO_proc_str'] = str(defining_process.shell_string(schannel=True,
                                        forbid=True, main=False, pdg_order=False, print_id = False) + '_')

        # write testR
        filename = pjoin(dirpath, 'testRR_%d_%d_%d_%d.f' %(isec,jsec,c3p,d3p) )
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

        # write testR
        filename = pjoin(dirpath, 'testRR_%d_%d_%d_%d.f' %(isec,jsec,ksec,lsec) )
        file = open(pjoin(dirmadnklo,"tmp_fortran/tmp_files/testRR_template.f")).read()
        file = file % replace_dict
        writer(filename).writelines(file)

        return True


    #===========================================================================
    # write driver_isec_jsec for real subprocess directory
    #===========================================================================

    def write_driver_npt_template(self, writer, dirpath, dirmadnklo, i , isec, jsec, c3p, d3p, UBgraphs):

        replace_dict = {}
        replace_dict['isec'] = isec
        replace_dict['jsec'] = jsec
        replace_dict['c3p'] = c3p
        replace_dict['d3p'] = d3p
        replace_dict['UBgraphs'] = UBgraphs

        # write driver
        filename = pjoin(dirpath, 'driver_RR_%d_%d_%d_%d.f' % (isec,jsec,c3p,d3p))
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
        proc_str += """PROC_FILES= get_UnderlyingProc_PDGs.o matrix_%s.o""" % defining_process.shell_string(
            schannel=True, forbid=True, main=False, pdg_order=False, print_id = False)

        replace_dict['proc_str'] = proc_str

        for i in range(0,len(overall_sector_info)):
            isec = overall_sector_info[i]['isec']
            jsec = overall_sector_info[i]['jsec']
            c3p = overall_sector_info[i]['c3p']
            d3p = overall_sector_info[i]['d3p']
            #ksec = overall_sector_info[i]['ksec']
            #lsec = overall_sector_info[i]['lsec']
            replace_dict['isec'] = isec
            replace_dict['jsec'] = jsec
            replace_dict['c3p'] = c3p
            replace_dict['d3p'] = d3p
            #replace_dict['ksec'] = ksec
            #replace_dict['lsec'] = lsec

            files_str += 'FILES_%d_%d_%d_%d= ' % (isec, jsec, c3p, d3p)
            # Born ME
            if (overall_sector_info[i]['Born_str'] == '' and overall_sector_info[i]['alt_Born_str'] == ''):
                files_str += '$(USR_FILES) matrix_' + 'dummy' + '.o '
                string = 'dummy'
            elif (overall_sector_info[i]['Born_str'] != '' and overall_sector_info[i]['alt_Born_str'] == ''):
                files_str += '$(USR_FILES) matrix_' + overall_sector_info[i]['Born_str'] + '.o '
                string = overall_sector_info[i]['Born_str']
            else:
                files_str += '$(USR_FILES) matrix_' + overall_sector_info[i]['alt_Born_str'] + '.o '
                string = overall_sector_info[i]['Born_str']

            files_str += 'configs_%s.o ' % string
            files_str += 'props_%s.o ' % string
            files_str += 'decayBW_%s.o ' % string
            files_str += 'leshouche_%s.o ' % string

            # Real ME
            if (overall_sector_info[i]['alt_Real_str'] == ''):
                files_str += 'matrix_' + overall_sector_info[i]['Real_str'] + '.o '
            else:
                files_str += 'matrix_' + overall_sector_info[i]['alt_Real_str'] + '.o '

            files_str += 'driver_RR_%d_%d_%d_%d.o ' % (isec, jsec, c3p, d3p)
            files_str += 'NNLO_RRsub_%d_%d_%d_%d.o ' % (isec, jsec, c3p, d3p)
            files_str += 'NNLO_IR_limits_%d_%d_%d_%d.o ' % (isec, jsec, c3p, d3p)
            files_str += 'testRR_%d_%d_%d_%d.o ' % (isec, jsec, c3p, d3p)
            files_str += 'NNLO_K_%d_%d_%d_%d.o ' % (isec, jsec, c3p, d3p)
            files_str += '$(PROC_FILES) $(COMMON_FILES)\n'
            all_str += ' sector_%d_%d_%d_%d' % (isec, jsec, c3p, d3p)
            sector_str += """
sector_%d_%d_%d_%d_libs: libs sector_%d_%d_%d_%d

sector_%d_%d_%d_%d: $(FILES_%d_%d_%d_%d)
\t$(DEFAULT_F_COMPILER) $(patsubst %%,$(OBJ)/%%,$(FILES_%d_%d_%d_%d)) $(LIBS) $(LIBSC) -o $@
""" %(isec, jsec, c3p, d3p, isec, jsec, c3p, d3p, isec, jsec, c3p, d3p, isec, jsec, c3p, d3p, isec, jsec, c3p, d3p)

        object_str = """
%.o: %.f $(INCLUDE)
\t$(DEFAULT_F_COMPILER) -c $(FFLAGS) $(FDEBUG) -o $(OBJ)/$@ $<

#%.o: $(PATH_TO_COMMON_FILES)/%.f $(INCLUDE)
#\t$(DEFAULT_F_COMPILER) -c $(FFLAGS) $(FDEBUG) -o $(OBJ)/$@ $<

%.o: $(PATH_TO_USR_FILES)/%.f90 $(INCLUDE)
	$(DEFAULT_F_COMPILER) -c $(FFLAGS) $(FDEBUG) -o $(OBJ)/$@ $<

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
        file = open(pjoin(dirmadnklo,"tmp_fortran/tmp_files/makefile_npt_template")).read()
        file = file % replace_dict
        writer(filename).write(file)

        return True


