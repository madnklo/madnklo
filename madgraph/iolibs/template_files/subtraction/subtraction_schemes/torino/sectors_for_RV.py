#
# Functions to build necessary template files for RV subprocess directories
#
import copy
import sys

import commons.generic_sectors as generic_sectors
import madgraph.various.misc as misc
import madgraph.core.diagram_generation as diagram_generation
import madgraph.fks.fks_common as fks_common
import madgraph.integrator.vectors as vectors
import logging
from collections import defaultdict

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
#  Necessary template files for real-virtual subprocess directories
#
#      ? - all_sector_list.inc
#      ? - NNLO_KRV_template.f
#      ? - NNLO_RVsub_template.f
#      ? - driver_npo_RV_template.f
#      ? - testRV_template.f
#      ? - NNLO_IR_limits_template.f
#      ? - get_Born_PDGs.f
#      ? - makefile_npo_RV_template.f
#      ? - virtual_recoiler.inc (needed for checking recoiler consistency)
#
#      ? - links from Born to Real-Virtual subproc directories
#
#==================================================================================


class SectorGeneratorRV(sectors.SectorGenerator):

    def write_RV_templates(self, contrib_definition, defining_process, counterterms, integrated_counterterms):

        model = defining_process.get('model')
        initial_state_PDGs, final_state_PDGs = defining_process.get_cached_initial_final_pdgs()
        all_PDGs = initial_state_PDGs, final_state_PDGs

        leglist = defining_process.get('legs')

        print('INTO RV SECTOR')
        print(defining_process.shell_string())
        print(contrib_definition.get_shell_name())
        print(contrib_definition.process_definition.get('id'))

        all_sectors = []
        all_sector_legs = []
        all_sector_id_legs = []
        all_sector_recoilers = []

        pert_dict = fks_common.find_pert_particles_interactions(model)
        colorlist = [model['particle_dict'][l['id']]['color'] for l in leglist]

        # Generate sectors:
        # list of possible singular configurations: from i and j to possible ij
        # (just 2-particle sectors as at NLO)

        ####################################################
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

                #if i is not a gluon, then j must not be a final state gluon
                if i['id'] != 21 and j['id'] == 21 and j['state']:
                    continue

                #if both i and j are gluons, then keep just the case in which i (number) < j (number)
                if i['id'] == 21 and j['id'] == 21 and j['state']:
                    if j.get('number') < i.get('number') :
                        continue

                #if j and i are quarks and antiquark in the final state, let j be the quark
                #  this is needed in order to comply with the fct combine_ij inside fks_common
                if i['id'] == -j['id'] and j['state']:
                    if j['id'] < 0:
                        continue

                ijlist = fks_common.combine_ij(fks_common.to_fks_leg(i, model),
                                               fks_common.to_fks_leg(j, model),
                                               model, pert_dict)

                # if not ijlist:
                #     ijlist = fks_common.combine_ij(fks_common.to_fks_leg(j, model),
                #                                    fks_common.to_fks_leg(i, model),
                #                                    model, pert_dict)

                for ij in ijlist:
                    # copy the defining process, remove i and j
                    # and replace them by ij.
                    new_process = copy.copy(defining_process)
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
                        logger.info('NNLO RV sector found, legs %d, %d' % a_sector['sector'].leg_numbers)

        if not all_sectors:
            logger.critical('WARNING, no sectors found for %s' % defining_process.nice_string())

        # up to here we have identified all the sectors.
        #  Now for each sector we need to find the corresponding counterterms
        #  We also need to add the variable all_sector_list to each of
        #  them, containing the list of all sectors (needed to normalise
        #  the weight function
        all_sector_list = [s['sector'].leg_numbers for s in all_sectors]
        all_sector_mass_list = [s['sector'].masses for s in all_sectors]
        all_sector_id_list = [s['sector'].id for s in all_sectors]

        print('All sectors for RV NNLO : ' + str(all_sector_list))
        print('All sectors id for RV NNLO : ' + str(all_sector_id_list))

        all_local_counterterms_list = []
        necessary_ct_list = [] # [S_RV_gi, S_RV_gj, C_RV_ij, SC_RV_ij]
        necessary_ct = [0] * (7*len(all_sectors))
        i = 0
        for s in all_sectors:
            s['sector'].all_sector_list = all_sector_list
            s['sector'].all_sector_mass_list = all_sector_mass_list
            s['sector'].all_sector_id_list = all_sector_id_list
            print('s all_sector_id_list : ' + str(all_sector_id_list))

            if counterterms is not None:
                s['counterterms'] = []
                necessary_ct_list_one = [0]*4
                for i_ct, ct in enumerate(counterterms):
                    #print('i_ct + ct : ' + str(i_ct) + ' and ' + str(ct))
                    #print('ct : ' + str(ct))
                    #print('masses : ' + str(all_sector_mass_list))
                    current = ct.nodes[0].current
                    singular_structure = current.get('singular_structure').substructures[0]
                    all_legs = singular_structure.get_all_legs()
                    if singular_structure.name()=='S':
                        if all_legs[0].n == s['sector'].leg_numbers[0]: # should match to "i"
                            s['counterterms'].append(i_ct)
                            necessary_ct_list_one[0] =  'S_RV_g' # 1
                            necessary_ct[i] = ct

                        if s['sector'].id[0] == 21 and s['sector'].id[1] == 21:
                            s['counterterms'].append(i_ct)
                            necessary_ct_list_one[1] = 'S_RV_g' # 1
                            necessary_ct[i+1] = ct

                    if singular_structure.name()=='C':
                        if s['sector'].all_sector_mass_list[1][-1] != 'ZERO':
                            print(str(s['sector'].all_sector_mass_list))
                            print(s['sector'].all_sector_mass_list[1][-1])
                            break

                        if not singular_structure.substructures:
                            # pure-collinear CT: include if the legs match those of the sector
                            if sorted([l.n for l in all_legs]) == sorted(s['sector'].leg_numbers):
                                s['counterterms'].append(i_ct)
                                if s['sector'].id[0] == 21 and s['sector'].id[1] == 21:
                                    necessary_ct_list_one[2] = 'HC_RV_gg' # 1
                                    #necessary_ct_list_one[3] = 'SC_RV_gg' # 1
                                elif s['sector'].id[0] == 21 and s['sector'].id[1] != 21:
                                    necessary_ct_list_one[2] = 'HC_RV_gq' # 1
                                    #necessary_ct_list_one[3] = 'SC_RV_gq' # 1
                                elif s['sector'].id[1] == 21 and s['sector'].id[0] != 21:
                                    necessary_ct_list_one[2] = 'HC_RV_gq' # 1
                                else :
                                    necessary_ct_list_one[2] = 'HC_RV_qqx' # 1
                                necessary_ct[i+2] = ct

                all_local_counterterms_list.append(s['counterterms'])
                necessary_ct_list.append(necessary_ct_list_one)

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

        print('necessary cts : ' + str(necessary_ct_list))


        # Find the corresponding integrated counterterms (Needed for RV)


######################################### Write fortran template files for n+1 body #############################################

        # Set writer
        writer = writers.FortranWriter

        dirmadnklo=os.getcwd()
        dirpath = pjoin(dirmadnklo,glob.glob("%s/NNLO_RV_x_R_*" % interface.user_dir_name[0])[0])
        # defining_process.shell_string() is not enough for loop processes "P1_#_xxxx"
        # solution: use dirpath to determine the loop index #
        # TODO: generalise if P1_1_xxxx and P1_2_xxxx exist (different loop index, same final state)
        subproc_base = pjoin(dirpath, 'SubProcesses')
        proc_string = defining_process.shell_string().split("_", 1)[1]
        dirpath = glob.glob(pjoin(subproc_base, "P1_*_%s" % proc_string))[0]

######### Import Born-level PDGs from proc/SupProcesses directory

        sys.path.append(pjoin(dirmadnklo,"%s/SubProcesses" % interface.user_dir_name[0]))
        import Born_PDGs as PDGs_from_Born
        import Real_PDGs as PDGs_from_Real

#==================================================================================
#  Necessary template files for real-virtual subprocess directories
#
#       - all_sector_list.inc
#       - NNLO_KRV_template.f
#       - NNLO_RVsub_template.f
#       - driver_npo_RV_template.f
#       - testRV_template.f
#       - NNLO_RV_IR_limits_template.f
#       - get_Born_PDGs.f
#       - makefile_npo_RV_template.f
#       - virtual_recoiler.inc (needed for checking recoiler consistency)
#
#       - links from Born to Real-Virtual subproc directories
#
#==================================================================================

######### Write all_sector_list.inc

        self.write_all_sector_list_include(writers.FortranWriter, dirpath, all_sector_list)
        len_sector_list = len(all_sector_list)
        K_sector_lists = defaultdict(lambda: defaultdict(list))

######### Write NNLO_KRV_isec_jsec.f, NNLO_RVsub_isec_jsec.f

        overall_sector_info = []
        # Set replace_dict for NNLO_KRV_isec_jsec.f
        replace_dict_ct = {}
        # Set replace_dict for NNLO_RVsub_isec_jsec.f
        replace_dict_int_real = {}
        replace_dict_limits = {}

        # ------------------------------------------------------------------------------
        # Initialize all counterterm placeholders to avoid KeyErrors in templates
        # This ensures even SC_gg / SC_gq (soft-collinear) terms have dummy entries.
        # They will be overwritten later when actual process names are found.
        # ------------------------------------------------------------------------------
        necessary_default_ct_list = ['S_RV_g', 'HC_RV_gg', 'HC_RV_gq', 'HC_RV_qqx']
        for ct in necessary_default_ct_list:
            replace_dict_limits['proc_prefix_%s' % ct] = 'dummy'
        # ------------------------------------------------------------------------------

        # List of necessary underlying Born and Virtual strings and particle PDGs
        Born_processes = []
        # List of dirpathLO of the necessary underlying Born
        path_Born_processes = []
        # Link LO files to each real process directory
        dirpathLO_head = pjoin(dirmadnklo,glob.glob("%s/LO_*" % interface.user_dir_name[0])[0])
        dirpathNLO_head = pjoin(dirmadnklo,glob.glob("%s/NLO_R_*" % interface.user_dir_name[0])[0])
        necessary_default_ct_list = ['S_RV_g', 'HC_RV_gg', 'HC_RV_gq', 'HC_RV_qqx']

        for i in range(0,len(all_sector_list)):
            list_M2 = []
            list_str_def_M2 = []
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
                'alt_Born_path' :   '',
                'Real_str'      :   '',
                'Virt_str'      :   '',
                'path_to_Real'  :   '',
                'path_to_Virt'  :   '',
            }
            sector_info['isec'] = isec
            sector_info['jsec'] = jsec
            sector_info['iref'] = iref

            # default mapping for final-state collinear kernels (abc) == (ijr)
            mapping = [('isec', isec), ('jsec', jsec), ('iref', iref)]
            sector_info['mapping'] = [mapping[0][1], mapping[1][1], mapping[2][1]]

            #specify (abc) mapping choice
            mapping_str = """ \
                iU = %s
                iS = %s
                iB = %s
                iA = 1 ! default azimuth for NLO
            """ % (mapping[0][0], mapping[1][0], mapping[2][0])
            overall_sector_info.append(sector_info)

            # Initialise NNLO_RV_IR_limits.f for every sector [ij]
            string = "c Collection of relevant limits for sector [%d,%d]" %(isec,jsec)
            NNLO_RV_IR_limits_tmp_path = dirmadnklo + '/tmp_fortran/tmp_files/NNLO_limits/'
            os.system('echo ' + string + ' > ' + NNLO_RV_IR_limits_tmp_path + 'IR_tmp.f')

            for j in range(0, len(necessary_ct_list[i])):
                if necessary_ct_list[i][j] ==  0:
                    continue
                elif j == 0:
                    if id_isec != 21:
                        raise MadEvent7Error('%d is not a gluon!' % isec)
                    list_str_def_M2.append('DOUBLE PRECISION M2_%s(-2:0)\n' % necessary_ct_list[i][j])
                    list_M2.append('CALL SUB_M2_%s(isec,xs,xp,wgt,xj,xjB,nitRV,1d0,wgt_chan,ierr,M2_%s)\n'
                                       % (necessary_ct_list[i][j], necessary_ct_list[i][j]))
                    list_M2.append('K%s=K%s+M2_%s\n'
                                       % (necessary_ct_list[i][j].split("_")[0], necessary_ct_list[i][j].split("_")[0], necessary_ct_list[i][j]))
                    list_M2.append('if(ierr.eq.1)goto 999\n')
                    # Write ct template in NNLO_RV_IR_limits
                    os.system('cat ' + NNLO_RV_IR_limits_tmp_path + '/' + necessary_ct_list[i][j] + '.f >> ' + NNLO_RV_IR_limits_tmp_path + 'IR_tmp.f')
                    K_sector_lists['S'][(isec,)].append((isec,jsec))
                elif j == 1:
                    continue
                elif j == 2:
                    if (isec == iref) or (jsec == iref):
                        raise MadEvent7Error('Wrong recoiler %d,%d,%d!' % (isec,jsec,iref))
                    list_str_def_M2.append('DOUBLE PRECISION M2_%s(-2:0)\n' % necessary_ct_list[i][j])
                    list_M2.append('CALL SUB_M2_%s(isec,jsec,iref,xs,xp,xsb,xpb,wgt,xj,nitRV,1d0,wgt_chan,ierr,M2_%s)\n'
                                   % (necessary_ct_list[i][j], necessary_ct_list[i][j]))
                    list_M2.append('K%s=K%s+M2_%s\n'
                                   % (necessary_ct_list[i][j].split("_")[0], necessary_ct_list[i][j].split("_")[0], necessary_ct_list[i][j]))
                    list_M2.append('if(ierr.eq.1)goto 999\n')
                    # Write ct template in NNLO_RV_IR_limits
                    os.system('cat ' + NNLO_RV_IR_limits_tmp_path + '/' + necessary_ct_list[i][j] + '.f >> ' + NNLO_RV_IR_limits_tmp_path + 'IR_tmp.f')
                    K_sector_lists['C'][(isec,jsec)].append((isec,iref))

            # outside loop on necessary_ct_list
            str_def_M2 = " ".join(list_str_def_M2)
            replace_dict_ct['str_def_M2'] = str_def_M2
            str_M2 = " ".join(list_M2)
            str_int_real = " ".join(list_int_real)
            replace_dict_ct['str_M2'] = str_M2
            replace_dict_int_real['str_int_real'] = str_int_real
            replace_dict_int_real['NLO_process'] = str(defining_process.shell_string(schannel=True,
                                        forbid=True, main=False, pdg_order=False, print_id = False))
            replace_dict_int_real['NNLO_RV_process'] = str(defining_process.shell_string(schannel=True,
                                        forbid=True, main=False, pdg_order=False, print_id = False))
            replace_dict_int_real['mapping_str'] = mapping_str
            replace_dict_int_real['NLO_proc_str'] = str(defining_process.shell_string(schannel=True,
                                        forbid=True, main=False, pdg_order=False, print_id = False) + '_')
            replace_dict_int_real['NNLO_RV_proc_str'] = str(defining_process.shell_string(schannel=True,
                                        forbid=True, main=False, pdg_order=False, print_id = False) + '_')

            if necessary_ct_list[i][0] != 0 or necessary_ct_list[i][1] != 0:
                # list of proc str permutations 'epem_ddx' for template
                uB_proc = necessary_ct[i*5].current.shell_string_user(
                            schannel=True, forbid=True, main=False, pdg_order=False, print_id = False)
                R_proc = str(defining_process.shell_string(schannel=True,
                                        forbid=True, main=False, pdg_order=False, print_id = False))

                # list of proc str permutations '1_epem_ddx' for directory
                uB_proc_str_1 = necessary_ct[i*5].current.shell_string_user()
                for j in range(0,len(uB_proc)):
                    dirpathLO = pjoin(dirpathLO_head, 'SubProcesses', "P1_%s" % uB_proc[j])
                    dirpathNLOR = pjoin(dirpathNLO_head, 'SubProcesses', "P1_%s" % R_proc)
                    dirpathNLOV = pjoin(dirmadnklo,glob.glob("%s/NLO_V_x_B_*" % interface.user_dir_name[0])[0])
                    Vsubproc_base = pjoin(dirpathNLOV, 'SubProcesses')
                    V_proc = uB_proc[j]
                    dirpathNLOV = glob.glob(pjoin(Vsubproc_base, "P1_*_%s" % V_proc))[0]
                    if os.path.exists(dirpathLO):
                        replace_dict_int_real['strUB'] = uB_proc[j]
                        replace_dict_limits['proc_prefix_S_RV_g'] = uB_proc[j]
                        overall_sector_info[i]['Born_str'] = uB_proc[j]
                        overall_sector_info[i]['Real_str'] = R_proc
                        overall_sector_info[i]['Virt_str'] = V_proc
                        overall_sector_info[i]['path_to_Born'] = dirpathLO
                        overall_sector_info[i]['path_to_Real'] = dirpathNLOR
                        overall_sector_info[i]['path_to_Virt'] = dirpathNLOV
                        if uB_proc[j] not in Born_processes:
                            Born_processes.append(uB_proc[j])
                            path_Born_processes.append(dirpathLO)
                        break
                    if j == len(uB_proc) - 1:
                        extra_uB_proc = uB_proc[0]
                        replace_dict_int_real['strUB'] = extra_uB_proc
                        replace_dict_limits['proc_prefix_S_RV_g'] = extra_uB_proc
                        overall_sector_info[i]['Born_str'] = extra_uB_proc

            if necessary_ct_list[i][2] != 0 :
                # Loop over sectors with final state particles only
                if isec > 2 and jsec > 2:
                    tmp_proc = 'proc_prefix_%s' % necessary_ct_list[i][2]
                    uB_proc = necessary_ct[i*5+2].current.shell_string_user(
                                schannel=True, forbid=True, main=False, pdg_order=False, print_id = False)
                    uB_proc_str_1 = necessary_ct[i*5+2].current.shell_string_user()
                    R_proc = str(defining_process.shell_string(schannel=True,
                                        forbid=True, main=False, pdg_order=False, print_id = False))
                    flag = False
                    for j in range(0,len(uB_proc)):
                        dirpathLO = pjoin(dirpathLO_head, 'SubProcesses', "P1_%s" % uB_proc[j])
                        dirpathNLOR = pjoin(dirpathNLO_head, 'SubProcesses', "P1_%s" % R_proc)
                        dirpathNLOV = pjoin(dirmadnklo,glob.glob("%s/NLO_V_x_B_*" % interface.user_dir_name[0])[0])
                        Vsubproc_base = pjoin(dirpathNLOV, 'SubProcesses')
                        V_proc = uB_proc[j]
                        dirpathNLOV = glob.glob(pjoin(Vsubproc_base, "P1_*_%s" % V_proc))[0]
                        if os.path.exists(dirpathLO):
                            replace_dict_int_real['strUB'] = uB_proc[j]
                            replace_dict_limits[tmp_proc] = uB_proc[j]
                            overall_sector_info[i]['Born_str'] = uB_proc[j]
                            overall_sector_info[i]['Real_str'] = R_proc
                            overall_sector_info[i]['Virt_str'] = V_proc
                            overall_sector_info[i]['path_to_Born'] = dirpathLO
                            overall_sector_info[i]['path_to_Real'] = dirpathNLOR
                            overall_sector_info[i]['path_to_Virt'] = dirpathNLOV
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
                                    overall_sector_info[i]['alt_Born_path'] = pjoin(dirpathLO_head, 'SubProcesses', "P1_%s"
                                                                                    % "_".join(['1',overall_sector_info[i]['alt_Born_str']]))
                                    #print(overall_sector_info[i]['alt_Born_str'])
                                    #print(overall_sector_info[i]['alt_Born_path'])

                                    flag = True
                                    break
                            if(flag == True):
                                break

            if necessary_ct_list[i][3] != 0 : #== 1:
                # Loop over sectors with final state particles only
                if isec > 2 and jsec > 2:
                    tmp_proc = 'proc_prefix_%s' % necessary_ct_list[i][3]
                    uB_proc = necessary_ct[i*5+2].current.shell_string_user(
                                schannel=True, forbid=True, main=False, pdg_order=False, print_id = False)
                    uB_proc_str_1 = necessary_ct[i*5+2].current.shell_string_user()
                    flag = False
                    for j in range(0, len(uB_proc)):
                        dirpathLO = pjoin(dirpathLO_head, 'SubProcesses', "P1_%s" % uB_proc[j])
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
                                    overall_sector_info[i]['alt_Born_path'] = pjoin(dirpathLO_head, 'SubProcesses', "P1_%s"
                                                                                    % "_".join(['1',overall_sector_info[i]['alt_Born_str']]))

                                    flag = True
                                    break
                            if(flag == True):
                                break
            ###

            overall_sector_info[i]['Born_PDGs'] = getattr(PDGs_from_Born, "leg_PDGs_%s" % overall_sector_info[i]['Born_str'])
            # write NNLO_RV_IR_limits
            filename = pjoin(dirpath, 'NNLO_RV_IR_limits_%d_%d.f' % (isec, jsec))
            if os.path.lexists("%s/V_proc_prefix.txt" % dirpath):
                os.remove("%s/V_proc_prefix.txt" % dirpath)
            os.symlink( "%s/proc_prefix.txt" % overall_sector_info[i]['path_to_Virt'], "%s/V_proc_prefix.txt" % dirpath )
            proc_V_pref = open(pjoin(dirpath,"V_proc_prefix.txt")).read()
            replace_dict_limits['V_long_proc_prefix'] = proc_V_pref
            file = open(NNLO_RV_IR_limits_tmp_path + 'IR_tmp.f').read()
            file = file % replace_dict_limits
            writer(filename).writelines(file)
            os.system('rm ' + NNLO_RV_IR_limits_tmp_path + 'IR_tmp.f')

            ###

            replace_dict_int_real['isec'] = isec
            replace_dict_int_real['jsec'] = jsec
            replace_dict_int_real['iref'] = iref
            replace_dict_int_real['UBgraphs'] = overall_sector_info[i]['Born_str']
            proc_RV_pref = open(pjoin(dirpath,"proc_prefix.txt")).read()
            replace_dict_int_real['long_proc_prefix'] = proc_RV_pref
            filename_int_real = pjoin(dirpath, 'NNLO_RVsub_%d_%d.f' % (isec, jsec))
            file_int_real = open(pjoin(dirmadnklo,"tmp_fortran/tmp_files/NNLO_RVsub_template.f")).read()
            file_int_real = file_int_real % replace_dict_int_real
            writer(filename_int_real).writelines(file_int_real)
            UBgraphs = overall_sector_info[i]['Born_str']
            self.write_driver_npo_rv_template(writer, dirpath, dirmadnklo, i , isec, jsec, UBgraphs)

            # write NNLO_KRV
            filename = pjoin(dirpath, 'NNLO_KRV_%d_%d.f' % (isec, jsec))
            file = open(pjoin(dirmadnklo,"tmp_fortran/tmp_files/NNLO_KRV_template.f")).read()
            file = file % replace_dict_ct
            writer(filename).writelines(file)

            # write NNLO_I1
            filename = pjoin(dirpath, 'NNLO_I1_%d_%d.f' % (isec, jsec))
            file = open(pjoin(dirmadnklo,"tmp_fortran/tmp_files/NNLO_I1_RV_template.f")).read()
            file = file % replace_dict_int_real
            writer(filename).writelines(file)

             # write NNLO_I12
            filename = pjoin(dirpath, 'NNLO_I12_%d_%d.f' % (isec, jsec))
            file = open(pjoin(dirmadnklo,"tmp_fortran/tmp_files/NNLO_I12_RV_template.f")).read()
            file = file % replace_dict_limits
            writer(filename).writelines(file)

            # write testRV
            self.write_testRV_template_file(writer, dirpath, dirmadnklo, defining_process,
                                                    i, isec, jsec, necessary_ct_list, mapping_str,all_sector_mass_list[i])

        # check on overall_sector_info length
        if len(overall_sector_info) != len(all_sector_list):
            raise MadEvent7Error('WARNING, the list of sector-dictionary entries is not compatible with the total number of sectors!')


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

        self.write_get_UnderLying_PDGs_file(writer, dirpath, overall_sector_info)

######### Write makefile_npo_RV_template

        self.write_makefile_rv_file(writers.FileWriter, dirpath, dirmadnklo, defining_process, overall_sector_info)

######### Write all_K_sector_list

        self.write_all_K_sector_list(writer,dirpath,len_sector_list,K_sector_lists)

######### Write ajob_isec_jsec

#        self.write_ajob_npo_file(writers.FileWriter, dirpath, dirmadnklo, overall_sector_info)

######### Link Born files to each real process directory

        self.link_files_from_B_R_V_to_RV_dir(dirpath, Born_processes, path_Born_processes, overall_sector_info)


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

            #if not glob.glob('%s/virtual_recoilers.inc' % dirpath):
            #    os.symlink( '%s/virtual_recoilers.inc' % dirpath_virtual, '%s/virtual_recoilers.inc' % dirpath)

            # writing virtual_recoilers.inc
            v_rec = recoiler_function.get_virtual_recoiler(getattr(PDGs_from_Real, "leg_PDGs_%s" % overall_sector_info[i]['Real_str']))
            data_v_rec = str(v_rec).replace('[','').replace(']','').replace(' ','').replace('(','').replace(')','')
            file = """ \
              integer, parameter :: len_iref = %d
              integer iref(2,len_iref)
              data iref/%s/
            """ % (len(v_rec), data_v_rec)
            filename = pjoin(dirpath, 'virtual_recoilers.inc')
            writer(filename).writelines(file)

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
    # write K_sector_list file
    #===========================================================================

    def write_all_K_sector_list(self,writer,dirpath,len_sector_list,K_sector_lists):

        file = """ \
          integer, parameter :: len  = %d
          integer l
          """ % (len_sector_list)

        minl = 3  # start from final state
        for type, entries in K_sector_lists.items():
            maxl = max(max(key) for key in entries.keys())
            ndims = len(next(iter(entries.keys())))
            if ndims == 1:
                file += """ \
          integer %s_SECTOR_LIST(%d:%d,LEN,2)
          """ % (type,minl,maxl)
            elif ndims == 2:
                file += """ \
          integer %s_SECTOR_LIST(%d:%d,%d:%d,LEN,2)
          """ % (type,minl,maxl,minl,maxl)

        file += """ \
        """

        for type, entries in K_sector_lists.items():

            ndims = len(next(iter(entries.keys())))
            file += """
!         data %s \n""" % (type)

            for key, lists in sorted(entries.items()):

                n_zeros = len_sector_list - len(lists)
                lists_extended = lists + [(0,0)]*n_zeros

                for n, (a,b) in enumerate(lists_extended, 1):
                    if n > len(lists):
                        if ndims == 1:
                            i = key[0]
                            file += """ \
          DATA (%s_SECTOR_LIST(%d,%d:%d,L),L=1,2) /%d*0/ \n""" % (type,i,n,len_sector_list,2*n_zeros)
                            break
                        elif ndims == 2:
                            i,j = key
                            file += """ \
          DATA (%s_SECTOR_LIST(%d,%d,%d:%d,L),L=1,2) /%d*0/ \n""" % (type,i,j,n,len_sector_list,2*n_zeros)
                            break
                    else:
                        if ndims == 1:
                            i = key[0]
                            file += """ \
          DATA (%s_SECTOR_LIST(%d,%d,L),L=1,2) /%d,%d/ \n""" % (type,i,n,a,b)
                        elif ndims == 2:
                            i,j = key
                            file += """ \
          DATA (%s_SECTOR_LIST(%d,%d,%d,L),L=1,2) /%d,%d/ \n""" % (type,i,j,n,a,b)

        filename = pjoin(dirpath, 'all_K_sector_list.inc')
        writer(filename).writelines(file)

        return True

    #===========================================================================
    # write driver_isec_jsec for real subprocess directory
    #===========================================================================

    def write_driver_npo_rv_template(self, writer, dirpath, dirmadnklo, i , isec, jsec, UBgraphs):

        replace_dict = {}
        replace_dict['isec'] = isec
        replace_dict['jsec'] = jsec
        replace_dict['UBgraphs'] = UBgraphs

        # write driver
        filename = pjoin(dirpath, 'driver_%d_%d.f' % (isec, jsec))
        file = open(pjoin(dirmadnklo,"tmp_fortran/tmp_files/driver_npo_RV_template.f")).read()
        file = file % replace_dict
        writer(filename).writelines(file)

        return True


    #===========================================================================
    # write file for testing limits, 'testRV.f'
    #===========================================================================

    def write_testRV_template_file(self, writer, dirpath, dirmadnklo, defining_process,
                                        i, isec, jsec, necessary_ct_list, mapping_str,mass_list):

        replace_dict = {}
        replace_dict['isec'] = isec
        replace_dict['jsec'] = jsec
        proc_RV_pref = open(pjoin(dirpath,"proc_prefix.txt")).read()
        replace_dict['long_proc_prefix'] = proc_RV_pref

        limit_str = ''
        if necessary_ct_list[i][0] != 0 : #Si limit
            limit_str += """
c
c     soft limit"""

            if mass_list[-1] != 'ZERO':
                limit_str += """
      e=[0d0,1d0]"""
            else:
                limit_str += """
      e=[1d0,1d0]
      l=[0d0,0d0]
      call do_limit_RV_%d_%d('Si      ',e,l)
"""%(isec,jsec)
        if necessary_ct_list[i][1] != 0 : #Sj limit
            #TODO for future: massive recoiler to be implemented
            limit_str += """
c
c     soft limit
      e=[1d0,1d0]
      l=[1d0,0d0]
      call do_limit_RV_%d_%d('Sj      ',e,l)
"""%(isec,jsec)
        # Loop over sectors with final state particles only
        if isec > 2 and jsec > 2:
            if necessary_ct_list[i][2] != 0 : #Cij
                limit_str += """
c
c     collinear limit
        e=[0d0,1d0]
        l=[0d0,0d0]
      call do_limit_RV_%d_%d('Cij     ',e,l)
"""%(isec,jsec)
            limit_str += """
c
c     spurious collinear limit
        e=[1d0,0d0]
        l=[0d0,0d0]
      call do_limit_RV_%d_%d('Cir     ',e,l)
"""%(isec,jsec)
            limit_str += """
c
c     spurious collinear limit
        e=[1d0,0d0]
        l=[1d0,0d0]
      call do_limit_RV_%d_%d('Cjr     ',e,l)
"""%(isec,jsec)
        elif isec > 2 and jsec <= 2:
            limit_str += """Collinear limits still to be specified in sectors.py """
            raise MadEvent7Error('Collinear limits still to be specified in sectors.py. ')

        replace_dict['limit_str'] = limit_str
        replace_dict['NLO_proc_str'] = str(defining_process.shell_string(schannel=True,
                                        forbid=True, main=False, pdg_order=False, print_id = False) + '_')
        replace_dict['NNLO_RV_proc_str'] = str(defining_process.shell_string(schannel=True,
                                        forbid=True, main=False, pdg_order=False, print_id = False) + '_')
        replace_dict['mapping_str'] = mapping_str

        # write testRV
        filename = pjoin(dirpath, 'testRV_%d_%d.f' %(isec,jsec) )
        file = open(pjoin(dirmadnklo,"tmp_fortran/tmp_files/testRV_template.f")).read()
        file = file % replace_dict
        writer(filename).writelines(file)

        return True


    #===========================================================================
    # write 'get_Born_PDGs.f' to find labels/flavours of n-body kinematics
    #===========================================================================

    def write_get_UnderLying_PDGs_file(self,writer, dirpath, overall_sector_info):

        file = ''
        file += """ \
          subroutine get_UnderLying_PDGs(isec,jsec,ksec,lsec,npart,UnderLying_leg_PDGs)
          implicit none
          include 'nexternal.inc'
          integer isec, jsec, ksec, lsec
          integer npart
          integer Born_leg_PDGs(nexternal_UB)
          integer UnderLying_leg_PDGs(npart)
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
          if(npart .eq. nexternal_UB) then
          UnderLying_leg_PDGs = Born_leg_PDGs
          else
          write(*,*) 'get_UnderLying_PDGs: error'
          write(*,*) 'npart must be equal to nexternal_UB'
          write(*,*) 'npart, nexternal_UB = ', npart, nexternal_UB
          write(*,*) 'exit...'
          stop
          endif
          return
          end
          """

        filename = pjoin(dirpath, 'get_UnderlyingProc_PDGs.f')
        writer(filename).writelines(file)

        return True



        filename = pjoin(dirpath, 'get_UnderlyingProc_PDGs.f')
        writer(filename).writelines(file)

        return True

    #===========================================================================
    # write 'makefile' for real-virtual subprocesses
    #===========================================================================

    def write_makefile_rv_file(self, writer, dirpath, dirmadnklo, defining_process, overall_sector_info):

        replace_dict = {}
        proc_str = ''
        files_str = ''
        sector_str = ''
        all_str = 'all: libs'
        # the RV matel is originally "loop_matrix.f" -- it is renamed "matrix_<str_proc>.f" for consistency
        os.rename("%s/RV_loop_matrix.f" % dirpath, "%s/matrix_%s.f" % (dirpath, defining_process.shell_string(
            schannel=True, forbid=True, main=False, pdg_order=False, print_id = False)))
        proc_str += """PROC_FILES= get_UnderlyingProc_PDGs.o matrix_%s.o matrix_R_%s.o """ % (
            defining_process.shell_string(
                schannel=True, forbid=True, main=False, pdg_order=False, print_id = False),
            defining_process.shell_string(
                schannel=True, forbid=True, main=False, pdg_order=False, print_id = False))

        seen_born = set()
        for item in overall_sector_info:
            born_str = item['Born_str']
            if item.get('path_to_Born') and born_str not in seen_born:
                proc_str += ' matrix_' + born_str + '.o'
                seen_born.add(born_str)

        replace_dict['proc_str'] = proc_str

        for i in range(0,len(overall_sector_info)):
            isec = overall_sector_info[i]['isec']
            jsec = overall_sector_info[i]['jsec']
            replace_dict['isec'] = isec
            replace_dict['jsec'] = jsec
            files_str += 'FILES_%d_%d= ' % (isec, jsec)
            files_str += '$(USR_FILES) driver_%d_%d.o ' % (isec, jsec)
            files_str += 'NNLO_RVsub_%d_%d.o ' % (isec, jsec)
            files_str += 'NNLO_RV_IR_limits_%d_%d.o ' % (isec, jsec)
            if not glob.glob("%s/matrix_%s.f" % (dirpath, overall_sector_info[i]['Born_str'])):
                files_str += 'configs_%s.o ' % overall_sector_info[i]['Born_str']
                files_str += 'props_%s.o ' % overall_sector_info[i]['Born_str']
                files_str += 'decayBW_%s.o ' % overall_sector_info[i]['Born_str']
                files_str += 'leshouche_%s.o ' % overall_sector_info[i]['Born_str']

            files_str += 'testRV_%d_%d.o ' % (isec, jsec)
            files_str += 'NNLO_KRV_%d_%d.o ' % (isec, jsec)
            files_str += 'NNLO_I1_%d_%d.o ' % (isec, jsec)
            files_str += 'NNLO_I12_%d_%d.o $(PROCESS) $(PROC_FILES) $(COMMON_FILES)\n' % (isec, jsec)
            all_str += ' sector_%d_%d' % (isec, jsec)
            sector_str += """
sector_%d_%d_libs: libs sector_%d_%d

sector_%d_%d: $(FILES_%d_%d)
\t$(DEFAULT_F_COMPILER) $(patsubst %%,$(OBJ)/%%,$(FILES_%d_%d)) $(LIBS) $(LIBSC) $(LINKLIBS) -o $@
""" %(isec, jsec,isec, jsec,isec, jsec,isec, jsec,isec,jsec)

        object_str = """
%.o: %.f $(INCLUDE) | $(OBJ)/RV_polynomial.o
\t$(DEFAULT_F_COMPILER) -c $(FFLAGS) $(FDEBUG) -o $(OBJ)/$@ $<

%.o: $(PATH_TO_USR_FILES)/%.f $(INCLUDE)
\t$(DEFAULT_F_COMPILER) -c $(FFLAGS) $(FDEBUG) -o $(OBJ)/$@ $<

%.o: $(PATH_TO_USR_FILES)/%.f90 $(INCLUDE)
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
        file = open(pjoin(dirmadnklo,"tmp_fortran/tmp_files/makefile_npo_RV_template")).read()
        file = file % replace_dict
        writer(filename).write(file)

        return True


    #===========================================================================
    # function for linking files to Real-Virtual subprocess directory
    #===========================================================================

    def link_files_from_B_R_V_to_RV_dir(self, dirpath, Born_processes, path_Born_processes, overall_sector_info):

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
                if not glob.glob("%s/matrix_%s.f" % (dirpath, overall_sector_info[i]['alt_Born_str'])):
                    os.symlink( "%s/matrix_%s.f" % (overall_sector_info[i]['alt_Born_path'], overall_sector_info[i]['alt_Born_str']),
                            "%s/matrix_%s.f" % (dirpath, overall_sector_info[i]['alt_Born_str']) )
                    os.symlink( overall_sector_info[i]['alt_Born_path'] + '/%s_spin_correlations.inc' % overall_sector_info[i]['alt_Born_str'],
                            dirpath + '/%s_spin_correlations.inc' % overall_sector_info[i]['alt_Born_str'] )
                continue

            if not glob.glob("%s/matrix_%s.f" % (dirpath, overall_sector_info[i]['Born_str'])):
                os.symlink( "%s/matrix_%s.f" % (overall_sector_info[i]['path_to_Born'], overall_sector_info[i]['Born_str']),
                            "%s/matrix_%s.f" % (dirpath, overall_sector_info[i]['Born_str']) )
                os.symlink( "%s/matrix_%s.f" % (overall_sector_info[i]['path_to_Real'], overall_sector_info[i]['Real_str']),
                            "%s/matrix_R_%s.f" % (dirpath, overall_sector_info[i]['Real_str']) )
                os.symlink( overall_sector_info[i]['path_to_Born'] + '/%s_spin_correlations.inc' % overall_sector_info[i]['Born_str'],
                            dirpath + '/%s_spin_correlations.inc' % overall_sector_info[i]['Born_str'] )
                os.symlink( overall_sector_info[i]['path_to_Real'] + '/%s_spin_correlations.inc' % overall_sector_info[i]['Real_str'],
                            dirpath + '/%s_spin_correlations.inc' % overall_sector_info[i]['Real_str'] )

                Vpref = open(pjoin(dirpath,"V_proc_prefix.txt")).read()
                os.symlink( "%s/V_polynomial.f" % overall_sector_info[i]['path_to_Virt'], "%s/V_polynomial.f" % dirpath )
                os.symlink( "%s/V_loop_matrix.f" % overall_sector_info[i]['path_to_Virt'], "%s/V_loop_matrix.f" % dirpath )
                os.symlink( "%s/V_improve_ps.f" % overall_sector_info[i]['path_to_Virt'], "%s/V_improve_ps.f" % dirpath )
                os.symlink( "%s/born_matrix.f" % overall_sector_info[i]['path_to_Virt'], "%s/V_born_matrix.f" % dirpath )
                os.symlink( "%s/V_CT_interface.f" % overall_sector_info[i]['path_to_Virt'], "%s/V_CT_interface.f" % dirpath )
                os.symlink( "%s/V_loop_num.f" % overall_sector_info[i]['path_to_Virt'], "%s/V_loop_num.f" % dirpath )
                os.symlink( "%s/helas_calls_ampb_1.f" % overall_sector_info[i]['path_to_Virt'], "%s/V_helas_calls_ampb_1.f" % dirpath )
                os.symlink( "%s/V_mp_compute_loop_coefs.f" % overall_sector_info[i]['path_to_Virt'], "%s/V_mp_compute_loop_coefs.f" % dirpath )
                os.symlink( "%s/mp_helas_calls_ampb_1.f" % overall_sector_info[i]['path_to_Virt'], "%s/V_mp_helas_calls_ampb_1.f" % dirpath )
                os.symlink( "%s/coef_construction_1.f" % overall_sector_info[i]['path_to_Virt'], "%s/V_coef_construction_1.f" % dirpath )
                os.symlink( "%s/loop_CT_calls_1.f" % overall_sector_info[i]['path_to_Virt'], "%s/V_loop_CT_calls_1.f" % dirpath )
                os.symlink( "%s/mp_coef_construction_1.f" % overall_sector_info[i]['path_to_Virt'], "%s/V_mp_coef_construction_1.f" % dirpath )
                os.symlink( "%s/V_TIR_interface.f" % overall_sector_info[i]['path_to_Virt'], "%s/V_TIR_interface.f" % dirpath )
                os.symlink( "%s/V_compute_color_flows.f" % overall_sector_info[i]['path_to_Virt'], "%s/V_compute_color_flows.f" % dirpath )
                os.symlink( "%s/V_user_access_subroutines.f" % overall_sector_info[i]['path_to_Virt'], "%s/V_user_access_subroutines.f" % dirpath )
                os.symlink( "%s/V_tir_cache_size.inc" % overall_sector_info[i]['path_to_Virt'], "%s/V_tir_cache_size.inc" % dirpath )
                os.symlink( "%s/V_process_info.inc" % overall_sector_info[i]['path_to_Virt'], "%s/V_process_info.inc" % dirpath )
                os.symlink( "%s/%sspin_correlations.inc" % (overall_sector_info[i]['path_to_Virt'], Vpref),"%s/%sspin_correlations.inc" % (dirpath, Vpref) )
                os.symlink( "%s/nsquaredSO.inc" % overall_sector_info[i]['path_to_Virt'], "%s/V_nsquaredSO.inc" % dirpath )
                os.symlink( "%s/ngraphs.inc" % overall_sector_info[i]['path_to_Virt'], "%s/V_ngraphs.inc" % dirpath )
                src_dir = "%s/MadLoop5_resources/" % overall_sector_info[i]['path_to_Virt']
                tgt_dir = "%s/MadLoop5_resources/" % dirpath
                for name in os.listdir(src_dir):
                    src_path = os.path.join(src_dir,name)
                    tgt_path = os.path.join(tgt_dir,name)
                    try:
                        os.remove(tgt_path)
                    except OSError:
                        pass
                    os.symlink(src_path, tgt_path)

        return #all_sectors
