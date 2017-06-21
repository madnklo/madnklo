################################################################################
#
# Copyright (c) 2009 The MadGraph5_aMC@NLO Development team and Contributors
#
# This file is a part of the MadGraph5_aMC@NLO project, an application which 
# automatically generates Feynman diagrams and matrix elements for arbitrary
# high-energy processes in the Standard Model and beyond.
#
# It is subject to the MadGraph5_aMC@NLO license which should accompany this 
# distribution.
#
# For more information, visit madgraph.phys.ucl.ac.be and amcatnlo.web.cern.ch
#
################################################################################

"""Methods and classes to export a list of contributions in the ME7_format."""

import glob
import logging
import os
import re
import shutil

import madgraph.core.base_objects as base_objects
import madgraph.core.color_algebra as color
import madgraph.core.helas_objects as helas_objects
import madgraph.core.contributions as contributions
import aloha as aloha
import madgraph.iolibs.files as files
import madgraph.iolibs.file_writers as writers
import madgraph.iolibs.template_files as template_files
import madgraph.various.banner as banner_mod
from madgraph import MadGraph5Error, InvalidCmd, MG5DIR
from madgraph.iolibs.files import cp, ln, mv

import madgraph.various.misc as misc

_file_path = os.path.split(os.path.dirname(os.path.realpath(__file__)))[0] + '/'
logger = logging.getLogger('madgraph.ME7Exporter')
pjoin = os.path.join

#===============================================================================
# ME7Exporter
#===============================================================================
class ME7Exporter(object):
    """Class to take care of exporting a set Contributions in the ME7 format.
    This class does not inherit from VirtualExporter because it is not responsible
    for exporting one set of amplitudes, but rather a list of Contributions."""

    def __init__(self, cmd_interface, noclean, group_subprocesses):
        """Initialize an ME7Exporter with an output path and a list of contributions"""
        
        self.group_subprocesses = group_subprocesses
        self.contributions      = cmd_interface._curr_contribs 
        self.export_dir         = cmd_interface._export_dir 
        self.options            = cmd_interface.options

        # Already initialize the exporter of each contribution
        for contrib in self.contributions:
            contrib.initialize_exporter(cmd_interface, noclean, group_subprocesses=group_subprocesses)

    def pass_information_from_cmd(self, cmd_interface):
        """ Pass information from the cmd_interface to this exporter and all contributions."""
        
        # For ME7, apply this to all individual contributions
        for contrib in self.contributions:
            contrib.pass_information_from_cmd(cmd_interface)

    def copy_template(self, model):
        """ Copy the template directories."""

        # Create the root directory if necessary
        if os.path.exists(self.export_dir):
            if noclean:
                raise InvalidCmd("The output path '%s' already exists. Clean it first.")
            else:
                shutil.rmtree(self.export_dir)
        os.makedirs(self.export_dir)

        # Forward the request for copying the template to each contribution
        for contrib in self.contributions:
            contrib.copy_template(model)
   
    def export(self, nojpeg, args=[]):
        """ Distribute and organize the export of all contributions. """
        
        # Forward the export request to each contribution
        for contrib in self.contributions:
            contrib.export(nojpeg=nojpeg, group_processes=self.group_subprocesses, args=args)

    def finalize(self, flaglist, interface_history):
        """ Distribute and organize the finalization of all contributions. """
        
        # Forward the finalize request to each contribution
        for contrib in self.contributions:
            # Must clean the aloha Kernel before each aloha export for each contribution
            aloha.aloha_lib.KERNEL.clean()
#            misc.sprint(contrib.nice_string())
#            if isinstance(contrib, contributions.Contribution_V):
#                misc.sprint( [[ p.nice_string() for p in me.get('processes')] for me in contrib.all_matrix_elements.get_matrix_elements()] )
            contrib.finalize(flaglist = flaglist, interface_history = interface_history)
