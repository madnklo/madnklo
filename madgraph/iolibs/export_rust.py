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
import madgraph.iolibs.file_writers as writers

import madgraph.various.misc as misc
from madgraph import MadGraph5Error, InvalidCmd, MG5DIR
from madgraph.iolibs.files import cp, ln, mv

import madgraph.various.misc as misc

_file_path = os.path.split(os.path.dirname(os.path.realpath(__file__)))[0]
logger = logging.getLogger('madgraph.RustExporter')
pjoin = os.path.join

class RustExporterException(MadGraph5Error):
    """ Daughter class specifying an exception occurring within the RustExporter. """
    pass

#===============================================================================
# RustExporter
#===============================================================================
class RustExporter(object):
    """Class to take care of exporting Rust source code for generating a low-level
     representation of MadNkLO integrands.
    This class does not inherit from VirtualExporter because it is not responsible
    for exporting one set of amplitudes, but rather a list of MadNkLO integrands."""

    static_template_path = pjoin(_file_path, os.pardir, 'Template', 'rust')
    dynamic_template_path = pjoin(_file_path, os.pardir, 'madgraph', 'iolibs', 'template_files', 'rust')

    def __init__(self, cmd_interface, export_options={}):
        """Initialize an Rust exporter with an instance of the main command line interface as well as additional
         export options which provide an access to all necessary information about how the rust backend should
         be exported, together with the list of MadNkLO integrands that must be exported."""

        self.cmd_interface      = cmd_interface
        self.export_dir         = cmd_interface._export_dir
        self.options            = cmd_interface.options
        self.model              = cmd_interface._curr_model
        self.export_options     = export_options

    def copy_template(self):
        """ Copy the template directories."""

        rust_export_path = pjoin(self.export_dir,'rust')

        # The overall exporter of MAdNkLO should have already done this copy of the
        # static rust template, so it does not need to be repeated here. But if more static
        # information needs to be added, then it can be done here.
        if not os.path.exists(rust_export_path):
            shutil.copytree(self.static_template_path, self.export_dir)
        # Create an empty integrands directory for hosting integrand implementation here
        # (an empty directory cannot be placed in the static_template_path because it cannot be tracked by git).
        os.makedirs(pjoin(rust_export_path,'integrands'))

    def export_global_resources(self, all_MEAccessors):
        """ Export global resources independent of any integrand or accessor, basically the base skeletton of the rust
        implementation. """

        # Copy the some static rust template files
        self.copy_template()

        # Below one should perform the export steps relating to rendering all entries of the
        # all_MEAccessors dictionary exposed to rust as well
        # TODO

        return

    def export(self, integrand):
        """ Export one particular integrand. """

        # The rust exporter currently only support LO integrands
        if integrand.contribution_definition.overall_correction_order.count('N')>0:
            return

        integrand_short_name = '%s_%d'%(integrand.get_short_name(),integrand.ID)
        integrand_export_path = pjoin(self.export_dir,'rust','integrands',integrand_short_name)
        # Create a directory specifically for this integrand
        os.makedirs(integrand_export_path)

        # Create in their the dynamically generate rust source code
        my_context = {
            'my_conditional_variable' : 4,
            'my_nested_conditional_boolean' : True,
            'myContextList' : ["A sheep", "A tulip", "PrintTHIS"],
        }
        my_replacement_dict = {
            'integrand_name' : integrand_short_name,
            'a_dynamical_statement' : 'println!("MadNkLO ftw!");'
        }
        rust_writer = writers.RustWriter(pjoin(integrand_export_path, 'integrand_%s.rs'%(integrand_short_name)),opt='w')
        rust_writer.writelines(
            open(pjoin(self.dynamic_template_path,'integrand_B.rs')).read(),
            context=my_context, replace_dictionary=my_replacement_dict
        )

    def compile(self):
        """ Compile the rust backend."""
        #TODO Note: this should be done as much as possible using makefile targets so that one can easily recompile
        # the entire distribution manually.
        return

    def finalize(self, all_MEAccessors, all_integrands):
        """Distribute and organize the finalization export of all accessors and integrands. """
        return
