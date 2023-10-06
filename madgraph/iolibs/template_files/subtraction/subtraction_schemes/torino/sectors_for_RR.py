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

    def write_RR_templates(self, model, initial_state_PDGs, final_state_PDGs, all_PDGs, leglist,
                            all_sectors, all_sector_legs, all_sector_id_legs, all_sector_recoiler,
                            all_sector_list, all_sector_mass_list, all_sector_id_list,
                            all_local_counterterms_list, necessary_ct_list, necessary_ct,
                            dirmadnklo, dirpath, defining_process):

######### Set writer
        writer = writers.FortranWriter

######### Write all_sector_list.inc
        self.write_all_sector_list_include(writers.FortranWriter, dirpath, all_sector_list)
    
        return all_sectors
