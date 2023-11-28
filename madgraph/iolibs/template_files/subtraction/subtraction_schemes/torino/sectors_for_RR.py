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
        #         [bq q q'], [bq q bq'] (0)
        #      4) [g g, g g],  (4)
        #         [g g, g q], [g g, g bq],  (3)
        #         [g g, bq q], [g q, g q], [g bq, g bq], [g q, g bq], [g q, g q'], [g bq, g bq'], [g q, g bq'],  (2)
        #         [g q, bq q], [g q, bq' q'],  (1)
        #         [bq q, bq q], [bq q, bq' q']  (0)

        # First divide according to the possibility of having 3 or 4 particle topologies

        ####################################################
        # 3-particle sectors
        fks_k_j_i = {}
        threep_sectors = []
        threep_sectors_id = []
        for i, col_i in zip(leglist, colorlist):
            if col_i == 1:
                continue
            if not i.get('state'):
                continue
            fks_k_j_i[i.get('number')] = [] # not strictly needed

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

                for k, col_k in zip(leglist, colorlist):
                    if col_k == 1:
                        continue
                    if k.get('number') == i.get('number') or k.get('number') == j.get('number'):
                        continue

        ####################################################


        # 4-particle sectors: NLO x NLO case
        fks_j_from_i = {}
        fks_l_from_k = {}
        fourp_sectors = []
        fourp_sectors_id = []
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
                        logger.info('First part of 4p NNLO sector found, legs %d, %d' % a_sector['sector'].leg_numbers)


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
                                #print('Leg number : ' + str(a_sector['sector']))
                                #gl
                                all_sector_legs.append(k.get('number'))
                                all_sector_legs.append(l.get('number'))
                                # keep track of the masses
                                a_sector['sector'].masses = (model.get('particle_dict')[k.get('id')]['mass'],
                                                     model.get('particle_dict')[l.get('id')]['mass'])
                                #print('Masses : ' + str(a_sector['sector'].masses))
# gl
                                # keep track of the particles' identity
                                a_sector['sector'].id = (k.get('id'), l.get('id'))
                                all_sector_id_legs.append(k.get('id'))
                                all_sector_id_legs.append(l.get('id'))
                                #print('Identities : ' + str(a_sector['sector'].id))

                                all_sectors.append(a_sector)
                                logger.info('Second part of 4p NNLO sector found, legs %d, %d' % a_sector['sector'].leg_numbers)

                        if len(kllist)==0:
                            continue 

                        # TODO: remove redundant sectors; in the symmetrised case
                        # Zijkl = ijkl + ijlk + jikl + jilk + klij + klji + lkij + lkji
                        tmp_sector = [i.get('number'),j.get('number'),k.get('number'),l.get('number')]
                        tmp_sector_id = [i.get('id'),j.get('id'),k.get('id'),l.get('id')]

                        if len(fourp_sectors) == 0 or tmp_sector not in combs:
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

        if not all_sectors:
            logger.critical('WARNING, no sectors found for %s' % defining_process.nice_string())
   
        print('All sectors for NNLO : ' + str(fourp_sectors))
        print('All sectors id for NNLO : ' + str(fourp_sectors_id))


        # Now for each sector we need to find the corresponding counterterms

        # Find the corresponding integrated counterterms (None for RR)



######### Set writer
        #writer = writers.FortranWriter

######### Write all_sector_list.inc
        #self.write_all_sector_list_include(writers.FortranWriter, dirpath, all_sector_list)
    
        return #all_sectors
