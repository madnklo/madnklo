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
import madgraph.interface.common_run_interface as common_run_interface
import aloha as aloha
import models.import_ufo as import_ufo
import models.model_reader as model_reader
import madgraph.iolibs.files as files
import madgraph.iolibs.file_writers as writers
import madgraph.iolibs.template_files as template_files
import madgraph.iolibs.export_v4 as export_v4
import madgraph.various.banner as banner_mod
from madgraph import MadGraph5Error, InvalidCmd, MG5DIR
from madgraph.iolibs.files import cp, ln, mv

import madgraph.various.misc as misc

_file_path = os.path.split(os.path.dirname(os.path.realpath(__file__)))[0] + '/'
logger = logging.getLogger('madgraph.ME7Exporter')
pjoin = os.path.join

template_path = pjoin(_file_path,os.pardir,'Template')

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
        self.model              = cmd_interface._curr_model

        # Already initialize the exporter of each contribution
        for contrib in self.contributions:
            contrib.initialize_exporter(cmd_interface, noclean, group_subprocesses=group_subprocesses)

    def pass_information_from_cmd(self, cmd_interface):
        """ Pass information from the cmd_interface to this exporter and all contributions."""
        
        # For ME7, apply this to all individual contributions
        for contrib in self.contributions:
            contrib.pass_information_from_cmd(cmd_interface)

    def update_make_opts(self):
        """ Synchronizes make_opts with the MG5/ME7 options."""
        make_opts = pjoin(self.export_dir, 'Source', 'make_opts')
    
        cpp_compiler = self.options['cpp_compiler'] if self.options['cpp_compiler'] else 'g++'
        fortran_compiler = self.options['fortran_compiler'] if self.options['fortran_compiler'] else 'gfortran'
        f2py_compiler = self.options['f2py_compiler'] if self.options['f2py_compiler'] else 'f2py'

        is_clang = misc.detect_if_cpp_compiler_is_clang(cpp_compiler)
        is_lc    = misc.detect_cpp_std_lib_dependence(cpp_compiler) == '-lc++'

        # list of the variable to set in the make_opts file
        for_update= {'MACFLAG' : '-mmacosx-version-min=10.7' if is_clang and is_lc else '',
                     'STDLIB' : '-lc++' if is_lc else '-lstdc++',
                     'STDLIB_FLAG' : '-stdlib=libc++' if is_lc and is_clang else '',
                     'DEFAULT_CPP_COMPILER' : cpp_compiler,
                     'DEFAULT_F_COMPILER' : fortran_compiler,
                     'DEFAULT_F2PY_COMPILER' : f2py_compiler}
    
        common_run_interface.CommonRunCmd.update_make_opts_full(make_opts, for_update)

    def copy_template(self, model):
        """ Copy the template directories."""

        # Create the root directory if necessary
        if os.path.exists(self.export_dir):
            if noclean:
                raise InvalidCmd("The output path '%s' already exists. Clean it first.")
            else:
                shutil.rmtree(self.export_dir)
        
        shutil.copytree(pjoin(template_path,'ME7'),self.export_dir)

        # Make sure to have make_opts synchronized
        self.update_make_opts()

        # Forward the request for copying the template to each contribution
        self.contributions.apply_method_to_all_contribs('copy_template', method_args = [model])
   
    def export(self, nojpeg, args=[]):
        """ Distribute and organize the export of all contributions. """
        
        # Forward the export request to each contribution
        self.contributions.apply_method_to_all_contribs('export', 
            method_args = [],
            method_opts = {'nojpeg':nojpeg, 'group_processes':self.group_subprocesses, 'args':args})

    def copy_model_resources(self):
        """Make the copy/symbolic links"""
        model_path = pjoin(self.export_dir, 'Source','MODEL')
        if os.path.exists(pjoin(model_path, 'ident_card.dat')):
            mv(pjoin(model_path,'ident_card.dat'), pjoin(self.export_dir, 'Cards'))
        cp(pjoin(model_path, 'param_card.dat'), pjoin(self.export_dir, 'Cards'))
        mv(pjoin(model_path, 'param_card.dat'), pjoin(self.export_dir, 'Cards','param_card_default.dat'))


    def dump_ME7(self, all_MEAccessors, all_integrands):
        """ Dumps all necessary information in order to bootstrap the ME7 interface.
        It is mostly all contained in all_MEAccessors and all_integrands."""
        
        open(pjoin(self.export_dir,"MadEvent7.db"),'w').write("TO BE FILLED")

    def create_run_card(self):
        """ Create the run card."""
        
        # For now use the default LO card
        run_card = banner_mod.RunCardME7()

        history = ''
        processes = [[v[0] for v in contrib.get_processes_map().values()] for contrib in self.contributions]
        proc_characteristic = {
            'ninitial':processes[0][0].get_ninitial(), 
            'loop_induced': len(self.contributions.get_contributions_of_type(contributions.Contribution_LIB)), 
            'colored_pdgs': range(1,7)+[21]}
        run_card.create_default_for_process(proc_characteristic, history, processes)

        run_card.write(pjoin(self.export_dir, 'Cards', 'run_card.dat'), 
            template=pjoin(self.export_dir, 'Cards', 'run_card.dat'), python_template=True )
        run_card.write(pjoin(self.export_dir, 'Cards', 'run_card_default.dat'), 
            template=pjoin(self.export_dir, 'Cards', 'run_card.dat'), python_template=True )

    def compile(self):
        """ Compile all contributions and the global ME7 resources (e.g. MODEL)"""
        
        # Compile the MODEL first
        if os.path.isdir(pjoin(self.export_dir,'Source','MODEL')):
            logger.info("Compiling global ME7 Model")
            misc.compile(arg=['../../lib/libmodel.a'], cwd=pjoin(self.export_dir,'Source','MODEL'), mode='fortran')
        
        # Compile all contributions
        self.contributions.apply_method_to_all_contribs('compile', log='Compiling')

    def finalize(self, flaglist, interface_history):
        """ Distribute and organize the finalization of all contributions. """
        
        # Make sure contributions are sorted at this stage
        self.contributions.sort_contributions()

        # Save all the global couplings to write out afterwards
        global_wanted_couplings = []
        # Forward the finalize request to each contribution
        for contrib in self.contributions:
            # Must clean the aloha Kernel before each aloha export for each contribution
            aloha.aloha_lib.KERNEL.clean()
#
#            misc.sprint(contrib.nice_string())
#            if isinstance(contrib, contributions.Contribution_V):
#                misc.sprint( [[ p.nice_string() for p in me.get('processes')] for me in 
#                                       contrib.all_matrix_elements.get_matrix_elements()] )
#
            wanted_couplings_to_add_to_global = contrib.finalize(flaglist = flaglist, 
                                                                 interface_history = interface_history)
            global_wanted_couplings.extend(wanted_couplings_to_add_to_global)

        # Generate the global ME7 MODEL
        if global_wanted_couplings:
            output_dir=pjoin(self.export_dir, 'Source', 'MODEL')
            # Writing out the model common to all the contributions that can share it
            model_export_options = { 'complex_mass'  : self.options['complex_mass_scheme'], 
                                     'export_format' : 'madloop', # So as to have access to lha_read_mp.f
                                     'mp'            : True, 
                                     'loop_induced'  : False }
            model_builder = export_v4.UFO_model_to_mg4(self.model, output_dir, model_export_options)
            model_builder.build(global_wanted_couplings)
        
        # Now possibly add content to the pool of global ME7 resources before removing superfluous files
        # and linking to necessary global ME7 resources
        for contrib in self.contributions:
            contrib.add_content_to_global_ME7_resources(self.export_dir)
            contrib.remove_superfluous_content()
            contrib.link_global_ME7_resources(self.export_dir)

        # Create the run_card
        self.create_run_card()

        # Add the cards generated in MODEL to the Cards directory
        self.copy_model_resources()
        # Now link the Sources files within each contribution
        for contrib in self.contributions:
            contrib.make_model_symbolic_link()
        
        # Copy the UFO model to the global ME7 resources Source directory
        ME7_ufo_path = pjoin(self.export_dir,'Sources','ME7_UFO_model_%s'%os.path.basename(self.model.get('modelpath')))
        shutil.copytree(self.model.get('modelpath'), ME7_ufo_path)
        # And clear compiled files in it
        for path in misc.glob(pjoin(ME7_ufo_path,'*.pkl'))+misc.glob(pjoin(ME7_ufo_path,'*.pyc')):
            os.remove(path)

        # Compile all contributions and global ME7 resources (e.g. MODEL)
        self.compile()
        
        # Now generate all the ME accessors and integrand.
        # Notice that some of the information provided here (RunCard, ModelReader, root_path, etc...)
        # can and will be overwritten by the actualized values when the ME7Interface will be launched.
        # We provide it here just so as to be complete.
        
        # Obtain all the accessors to the Matrix Element made available in this process output
        all_MEAccessors = contributions.MEAccessorDict()
        for contrib in self.contributions:
            all_MEAccessors.add_MEAccessors(contrib.get_MEAccessors(self.export_dir))
        
        # Now generate all the integrands from the contributions exported
        all_integrands = []
        run_card = banner_mod.RunCardME7(pjoin(self.export_dir,'Cards','run_card.dat'))

        # We might want to recover whether prefix was used when importing the model and whether
        # the MG5 name conventions was used. But this is a detail that can easily be fixed later.
        modelReader_instance = import_ufo.import_model(
            pjoin(self.export_dir,'Sources','ME7_UFO_model_')+self.model.get('name'), prefix=True,
            complex_mass_scheme=self.options['complex_mass_scheme'] )
        modelReader_instance.pass_particles_name_in_mg_default()
        modelReader_instance.set_parameters_and_couplings(
                param_card = pjoin(self.export_dir,'Cards','param_card.dat'), 
                scale=run_card['scale'], 
                complex_mass_scheme=self.options['complex_mass_scheme'])
        
        for contrib in self.contributions:
            all_integrands.extend(contrib.get_integrands( 
                                        modelReader_instance, run_card, all_MEAccessors, self.options))
        
        # And finally dump ME7 output information so that all relevant objects
        # can be reconstructed for a future launch with ME7Interface.
        # Normally all the relevant information should simply be encoded in only:
        #  'all_MEAccessors' and 'all_integrands'.
        self.dump_ME7(all_MEAccessors, all_integrands)
        
        #
        # WARNING THE COE BELOW IS JUST FOR TESTING PURPOSES
        #
        
        # This is now just for gigs. Integrate that beast!
        # Of course, what should really happen is that the users starts a ME7_interface, that 
        # bootstraps from the dump above and starts the integration below with lunch.
        # So this is really just for testing purposes.
        import madgraph.integrator.integrators as integrators
        integrator_naive = integrators.SimpleMonteCarloIntegrator(all_integrands,
            **{'n_iterations'            : 10,
               'n_points_per_iterations' : 100,
               'accuracy_target'         : None,
               'verbosity'               : 1 }
            )
        import madgraph.integrator.pyCubaIntegrator as pyCubaIntegrator        
        integrator_vegas = pyCubaIntegrator.pyCubaIntegrator(all_integrands, 
            **{'algorithm' : 'Vegas', 
               'verbosity' : 1,
               'seed'      : 3,
               'target_accuracy' : 1.0e-3,
               'n_start'   : 1000,
               'n_increase': 500,
               'n_batch'   : 1000,
               'max_eval'  : 100000,
               'min_eval'  : 0})

        # Now run them all!
        for integrator in [integrator_naive, integrator_vegas]:
            xsec, error = integrator.integrate()
            logger.info("="*100)
            logger.info('{:^100}'.format("\033[92mCross-section for process output '%s' with integrator '%s':\033[0m"
                                                                 %(self.export_dir, integrator.get_name())))
            logger.info('{:^100}'.format("\033[94m%.5e +/- %.2e [pb]\033[0m"%(xsec, error)))
            logger.info("="*100+"\n")
        