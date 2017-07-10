################################################################################
#
# Copyright (c) 2011 The MadGraph5_aMC@NLO Development team and Contributors
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
"""A user friendly command line interface to steer ME7 integration.
   Uses the cmd package for command interpretation and tab completion.
"""
from __future__ import division

import collections
import glob
import logging
import math
import os
import random
import re
import math

import stat
import subprocess
import sys
import time
import tarfile
import StringIO
import shutil
import copy

try:
    import readline
    GNU_SPLITTING = ('GNU' in readline.__doc__)
except:
    GNU_SPLITTING = True

root_path = os.path.split(os.path.dirname(os.path.realpath( __file__ )))[0]
root_path = os.path.split(root_path)[0]
sys.path.insert(0, os.path.join(root_path,'bin'))

# useful shortcut
pjoin = os.path.join
# Special logger for the Cmd Interface
logger = logging.getLogger('madevent7.stdout') # -> stdout
logger_stderr = logging.getLogger('madevent7.stderr') # ->stderr

try:
    import madgraph
except ImportError: 
    # import from madevent directory
    MADEVENT7 = True
    import internal.extended_cmd as cmd
    import internal.common_run_interface as common_run
    import internal.banner as banner_mod
    import internal.misc as misc
    from internal import InvalidCmd, MadGraph5Error, ReadWrite
    import internal.files as files
    import internal.save_load_object as save_load_object
    import internal.cluster as cluster
    import internal.check_param_card as check_param_card
    import internal.model_reader as model_reader
    import internal.import_ufo as import_ufo
    import internal.lhe_parser as lhe_parser
    import internal.integrands as integrands
    import internal.integrators as integrators
    import internal.phase_space_generators as phase_space_generators
    import internal.pyCubaIntegrator as pyCubaIntegrator
    import internal.vegas3_integrator as vegas3_integrator

#    import internal.histograms as histograms # imported later to not slow down the loading of the code
    from internal.files import ln
else:
    # import from madgraph directory
    MADEVENT7 = False
    import madgraph.interface.extended_cmd as cmd
    import madgraph.interface.common_run_interface as common_run
    import madgraph.iolibs.files as files
    import madgraph.iolibs.save_load_object as save_load_object
    import madgraph.various.banner as banner_mod
    import madgraph.various.cluster as cluster
    import madgraph.various.misc as misc
    import madgraph.various.lhe_parser as lhe_parser
    import madgraph.integrator.integrands as integrands
    import madgraph.integrator.integrators as integrators
    import madgraph.integrator.phase_space_generators as phase_space_generators
    import madgraph.integrator.pyCubaIntegrator as pyCubaIntegrator
    import madgraph.integrator.vegas3_integrator as vegas3_integrator
#    import madgraph.various.histograms as histograms  # imported later to not slow down the loading of the code
    import models.check_param_card as check_param_card
    import models.model_reader as model_reader
    import models.import_ufo as import_ufo

    from madgraph.iolibs.files import ln    
    from madgraph import InvalidCmd, MadGraph5Error, MG5DIR, ReadWrite



class MadEvent7Error(Exception):
    pass

class ZeroResult(MadEvent7Error):
    pass

#===============================================================================
# CmdExtended
#===============================================================================
class CmdExtended(common_run.CommonRunCmd):
    """Particularisation of the cmd command for MadEvent7"""

    #suggested list of command
    next_possibility = {
        'start': [],
    }
    
    debug_output = 'ME7_debug'
    error_debug = 'Please report this bug on https://bugs.launchpad.net/mg5amcnlo\n'
    error_debug += 'More information is found in \'%(debug)s\'.\n' 
    error_debug += 'Please attach this file to your report.'

    config_debug = 'If you need help with this issue please contact us on https://answers.launchpad.net/mg5amcnlo\n'


    keyboard_stop_msg = """stopping all operation
            in order to quit MadGraph5_aMC@NLO please enter exit"""
    
    # Define the Error
    InvalidCmd = InvalidCmd
    ConfigurationError = MadGraph5Error

    def __init__(self, me_dir, options, *arg, **opt):
        """Init history and line continuation"""
        
        # Tag allowing/forbiding question
        self.force = False
        
        # If possible, build an info line with current version number 
        # and date, from the VERSION text file
        info = misc.get_pkg_info()
        info_line = ""
        if info and info.has_key('version') and  info.has_key('date'):
            len_version = len(info['version'])
            len_date = len(info['date'])
            if len_version + len_date < 30:
                info_line = "#*         VERSION %s %s %s         *\n" % \
                            (info['version'],
                            (30 - len_version - len_date) * ' ',
                            info['date'])
        else:
            version = open(pjoin(root_path,'MGMEVersion.txt')).readline().strip()
            info_line = "#*         VERSION %s %s                *\n" % \
                            (version, (24 - len(version)) * ' ')    

        # Create a header for the history file.
        # Remember to fill in time at writeout time!
        self.history_header = \
        '#************************************************************\n' + \
        '#*               MadGraph5_aMC@NLO/MadEvent                 *\n' + \
        '#*                                                          *\n' + \
        "#*                *                       *                 *\n" + \
        "#*                  *        * *        *                   *\n" + \
        "#*                    * * * * 5 * * * *                     *\n" + \
        "#*                  *        * *        *                   *\n" + \
        "#*                *                       *                 *\n" + \
        "#*                                                          *\n" + \
        "#*                                                          *\n" + \
        info_line + \
        "#*                                                          *\n" + \
        "#*    The MadGraph5_aMC@NLO Development Team - Find us at   *\n" + \
        "#*    https://server06.fynu.ucl.ac.be/projects/madgraph     *\n" + \
        '#*                                                          *\n' + \
        '#************************************************************\n' + \
        '#*                                                          *\n' + \
        '#*               Command File for MadEvent                  *\n' + \
        '#*                                                          *\n' + \
        '#*     run as ./bin/madevent.py filename                    *\n' + \
        '#*                                                          *\n' + \
        '#************************************************************\n'
        
        if info_line:
            info_line = info_line[1:]

        logger.info(\
        "************************************************************\n" + \
        "*                                                          *\n" + \
        "*                     M A D E V E N T 7                    *\n" + \
        "*                                                          *\n" + \
        "*                 *                       *                *\n" + \
        "*                   *        * *        *                  *\n" + \
        "*                     * * * * 5 * * * *                    *\n" + \
        "*                   *        * *        *                  *\n" + \
        "*                 *                       *                *\n" + \
        "*                                                          *\n" + \
        info_line + \
        "*                                                          *\n" + \
        "*    The MadGraph5_aMC@NLO Development Team - Find us at   *\n" + \
        "*    https://server06.fynu.ucl.ac.be/projects/madgraph     *\n" + \
        "*                                                          *\n" + \
        "*               Type 'help' for in-line help.              *\n" + \
        "*                                                          *\n" + \
        "************************************************************")
        super(CmdExtended, self).__init__(me_dir, options, *arg, **opt)
        
    def get_history_header(self):
        """return the history header""" 
        return self.history_header % misc.get_time_info()        
        
#===============================================================================
# HelpToCmd
#===============================================================================
class HelpToCmd(object):
    """ The Series of help routine for MadEvent7Cmd"""
     
    def help_launch(self):
        """ TODO"""
        logger.info("syntax: launch [run_name] [options])")

#===============================================================================
# CheckValidForCmd
#===============================================================================
class ParseCmdArguments(object):
    """ The Series of check routine for MadEvent7Cmd"""
    
    def parse_launch(self, args):
        """ Parses the argument of the launch command."""

        launch_options = {'integrator': 'VEGAS3',
                          'n_points': None,
                          'n_iterations':None,
                          'verbosity':0}        
        #for name in ['Naive','VEGAS','VEGAS3','SUAVE','DIVONNE','CUHRE']:
        
        for arg in args:
            try:
                key, value = arg.split('=',1)
            except ValueError:
                key = arg
                value = None
            
            if key == '--integrator':
                if value not in self._integrators:
                    raise InvalidCmd("Selected integrator '%s' not reckognized."%value)
                launch_options['integrator'] = value
            elif key in ['--n_points', '--n_iterations']:
                launch_options[key[2:]] = int(value)
            elif key=='--verbosity':
                modes = {'none':0, 'integrator':1, 'all':2}
                launch_options[key[2:]] = modes[value.lower()]
            else:
                raise InvalidCmd("Option '%s' for the launch command not reckognized."%key)

        return launch_options

#===============================================================================
# CompleteForCmd
#===============================================================================
class CompleteForCmd(cmd.CompleteCmd):
    """ The Series of help routine for MadEvent7Cmd"""

    def complete_launch(self, *args, **opts):
        # TODO
        return []

#===============================================================================
# MadEvent7Cmd
#===============================================================================
class MadEvent7Cmd(CompleteForCmd, CmdExtended, ParseCmdArguments, HelpToCmd, common_run.CommonRunCmd):
    """The command line processor for MadEvent7"""    

    _set_options = []
    
    # For now a very trivial setup, but we will want of course to have all these meta
    # parameter controllable by user commands, eventually.
    integrator_verbosity = 0 if logger.level > logging.DEBUG else 1
    _integrators = {
       'NAIVE' : (integrators.SimpleMonteCarloIntegrator, 
                  { 'n_iterations'            : 10,
                    'n_points_per_iterations' : 100,
                    'accuracy_target'         : None,
                    'verbosity'               : integrator_verbosity  } ),
    
       'VEGAS3' : (vegas3_integrator.Vegas3Integrator,
                   { 'survey_n_iterations'     : 10,
                     'survey_n_points'         : 1000,
                     'refine_n_iterations'     : 10,
                     'refine_n_points'         : 2000,
                     'verbosity'               : integrator_verbosity  } ),
    
       'VEGAS' : (pyCubaIntegrator.pyCubaIntegrator, 
                  { 'algorithm' : 'Vegas', 
                    'verbosity' : integrator_verbosity,
                    'seed'      : 3,
                    'target_accuracy' : 1.0e-3,
                    'n_start'   : 1000,
                    'n_increase': 500,
                    'n_batch'   : 1000,
                    'max_eval'  : 100000,
                    'min_eval'  : 0}),
       
       'SUAVE'   : (pyCubaIntegrator.pyCubaIntegrator, 
                    { 'algorithm' :'Suave', 
                      'verbosity' : integrator_verbosity,
                      'target_accuracy' : 1.0e-3,
                      'max_eval'  : 100000,
                      'min_eval'  : 0 } ),
      
       'DIVONNE' : (pyCubaIntegrator.pyCubaIntegrator, 
                    { 'algorithm' : 'Divonne', 
                      'verbosity': integrator_verbosity,
                      'target_accuracy' : 1.0e-5,
                      'max_eval'  : 100000000,
                      'min_eval'  : 0 } ),
    
       'CUHRE'   : (pyCubaIntegrator.pyCubaIntegrator, 
                    { 'algorithm' : 'Cuhre',
                      'verbosity' : integrator_verbosity,
                      'target_accuracy' : 1.0e-3,
                      'max_eval'  : 100000,
                      'min_eval'  : 0 } ),
    }
    
    def __init__(self, me_dir = None, options={}, *completekey, **stdin):
        """ Initialize the interface """

        # Temporary n_initial value which will be updated during the bootstrap
        self.n_initial = 2
        CmdExtended.__init__(self, me_dir, options, *completekey, **stdin)
        self.prompt = "ME7 @ %s > "%os.path.basename(pjoin(self.me_dir))
        
        # Initialize default properties that will be overwritten during the bootstrap.
        self.all_MEAccessors = None
        self.all_integrands = []
        self.model = None
        self.run_card = None
        # Overwrite the above properties from the bootstrap using the database written at the output stage
        self.bootstrap_ME7()
        
        # Instance of the current integrator
        self.integrator = None
        
    def check_output_type(self, path):
        """ Check that the output path is a valid Madevent 7directory """        
        return os.path.isfile(os.path.join(path,'MadEvent7.db'))

    def check_already_running(sefl):
        pass

    def set_configuration(self, amcatnlo=False, final=True, **opt):
        """ Assign all configuration variable from various files."""
        
        # This will basically only load the options
        super(MadEvent7Cmd,self).set_configuration(amcatnlo=False, final=final, **opt)

    def bootstrap_ME7(self):
        """ Setup internal state of MadEvent7, mostly the integrand, from the content of
        what was dumped in the database MadEvent7.db at the output stage."""
        
        # Load the current run_card
        self.run_card = banner_mod.RunCardME7(pjoin(self.me_dir,'Cards','run_card.dat'))
        
        # Now reconstruct all the relevant information from ME7_dump
        logger.info("Loading MadEvent7 setup from file '%s'."%pjoin(self.me_dir,'MadEvent7.db'))
        ME7_dump = save_load_object.load_from_file(pjoin(self.me_dir,'MadEvent7.db'))
                
        # We might want to recover whether prefix was used when importing the model and whether
        # the MG5 name conventions was used. But this is a detail that can easily be fixed later.
        self.model = import_ufo.import_model(
            pjoin(self.me_dir,'Source',ME7_dump['model_name']), prefix=True,
            complex_mass_scheme = ME7_dump['model_with_CMS'] )
        self.model.pass_particles_name_in_mg_default()
        self.model.set_parameters_and_couplings(
                param_card = pjoin(self.me_dir,'Cards','param_card.dat'), 
                scale=self.run_card['scale'], 
                complex_mass_scheme=ME7_dump['model_with_CMS'])

        self.all_MEAccessors = ME7_dump['all_MEAccessors']['class'].initialize_from_dump(
                                                ME7_dump['all_MEAccessors'], root_path = self.me_dir)
        self.all_integrands = [integrand_dump['class'].initialize_from_dump(
            integrand_dump, self.model, self.run_card, self.all_MEAccessors, self.options
                                                    ) for integrand_dump in ME7_dump['all_integrands']]

        self.n_initial = ME7_dump['n_initial']
        
    def compile(self):
        """ Re-compile all necessary resources, so as to make sure it is up to date."""
        # TODO
        pass

    def do_launch(self, line, *args, **opt):
        """Main command, starts the cross-section computation. Very basic setup for now.
        We will eventually want to have all of these meta-data controllable via user commands
        in the interface.
        This is super naive and only for illustrative purposes for now."""
        
        args = self.split_arg(line)
        
        launch_options = self.parse_launch(args)

        # In principle we want to start by recompiling the process output so as to make sure
        # that everything is up to date.
        self.compile()

        # Re-initialize the integrator
        integrator_name = launch_options['integrator']
        integrator_options = self._integrators[integrator_name][1]
        
        integrator_options['verbosity'] = launch_options['verbosity']
        
        if launch_options['n_points']:
            if integrator_name=='VEGAS3':
                integrator_options['survey_n_points'] = launch_options['n_points']
                integrator_options['refine_n_points'] = 5*launch_options['n_points']
            elif integrator_name=='NAIVE':
                integrator_options['n_points_per_iterations'] = launch_options['n_points']
                
            else:
                # For now support this option only for some integrators
                raise InvalidCmd("The options 'n_points' is not supported for the integrator %s."%integrator_name)

        if launch_options['n_iterations']:
            if integrator_name=='VEGAS3':
                integrator_options['survey_n_iterations'] = launch_options['n_iterations']
                integrator_options['refine_n_iterations'] = launch_options['n_iterations']
            elif integrator_name=='NAIVE':
                integrator_options['n_iterations'] = launch_options['n_iterations']
            else:
                # For now support this option only for some integrators
                raise InvalidCmd("The options 'n_iterations' is not supported for the integrator %s."%integrator_name)

        self.integrator = self._integrators[integrator_name][0](self.all_integrands, **integrator_options)
        
        if len(set([len(itgd.get_dimensions()) for itgd in self.all_integrands]))>1 and integrator_name not in ['NAIVE']:
            # Skip integrators that do not support integrands with different dimensions.
            # Of course, in the end we wil not naively automatically put all integrands alltogether in the same integrator
            raise InvalidCmd("For now, whenever you have several integrands of different dimensions, only "+
                             "the NAIVE integrator is available.")
        
        # The integration can be quite verbose, so temporarily setting their level to 50 by default is best here
        if launch_options['verbosity'] > 1:
            logger_level = logging.DEBUG
        else:
            logger_level = logging.INFO
        
        with misc.MuteLogger(['madevent7.stdout','madevent7.stderr'],[logger_level,logger_level]):
            xsec, error = self.integrator.integrate()
            
        logger.info("="*100)
        logger.info('{:^100}'.format("\033[92mCross-section with integrator '%s':\033[0m"%self.integrator.get_name()))
        logger.info('{:^100}'.format("\033[94m%.5e +/- %.2e [pb]\033[0m"%(xsec, error)))
        logger.info("="*100+"\n")
        
    def do_show_grid(self, line):
        """ Minimal implementation for now with no possibility of passing options."""
        
        show_grid_options = {'n_grid':40,
                             'shrink':False,
                             'axes':None}
        
        if hasattr(self.integrator, 'show_grid'):
            self.integrator.show_grid(**show_grid_options)
        else:
            raise InvalidCmd('The current integrator used (%s) cannot display grids.'%self.integrator.get_name())
        
    def get_characteristics(self, path=None):
        """reads process characteristics and initializes the corresponding dictionary"""

        # ME7 interface doesn't use this structure much, but some minimal information is
        # still necessary for some of the parent's class features to work OK.
        self.proc_characteristics = {'ninitial': self.n_initial}
        return self.proc_characteristics
  
  
    def configure_directory(self):
        """ All action require before any type of run """   
        misc.sprint("Configure directory -- Nothing to be done here yet.")        
        pass

#===============================================================================
# MadEvent7CmdShell
#===============================================================================
class MadEvent7CmdShell(MadEvent7Cmd, cmd.CmdShell):
    """The shell command line processor of MadEvent7"""  
    
    
#===============================================================================
# ME7Integrand
#===============================================================================    
class ME7Integrand(integrands.VirtualIntegrand):
    """ Specialization for multi-purpose integration with ME7."""
    
    # This parameter defines a thin layer around the boundary of the unit hypercube of the 
    # random variables generating the phase-space, so as to avoid extrema which are an issue in most
    # PS generators.
    epsilon_border = 1e-10
    
    # Maximum size of the cache for PDF calls
    PDF_cache_max_size = 1000
    
    # The lowest value that the center of mass energy can take.
    # We take here 1 GeV, as anyway below this non-perturbative effects dominate and factorization does not
    # make sense anymore
    absolute_Ecm_min = 1.
    
    def __new__(cls, model, 
                     run_card,
                     contribution_definition,
                     processes_map,
                     all_MEAccessors,
                     ME7_configuration, **opt):
        all_args = [model, run_card, contribution_definition, processes_map,
                    all_MEAccessors, ME7_configuration ]
        if cls is ME7Integrand:
            target_class = None
            if contribution_definition.correction_order == 'LO':
                if contribution_definition.n_loops == 0 and \
                  contribution_definition.n_unresolved_particles == 0:
                   target_class = ME7Integrand_B
                elif contribution_definition.n_loops == 1 and \
                  contribution_definition.n_unresolved_particles == 0:
                   target_class = ME7Integrand_LIB
            elif contribution_definition.correction_order == 'NLO':
                if contribution_definition.n_loops == 1 and \
                   contribution_definition.n_unresolved_particles == 0:
                    target_class = ME7Integrand_V
                elif contribution_definition.n_loops == 0 and \
                     contribution_definition.n_unresolved_particles == 1:
                    target_class = ME7Integrand_R
            elif contribution_definition.correction_order == 'NNLO':
                if contribution_definition.n_loops == 0 and \
                   contribution_definition.n_unresolved_particles == 2:
                    target_class = ME7Integrand_RR                    
                else:
                    raise MadGraph5Error("Some NNLO type of integrands are not implemented yet.")                
            if not target_class:
                raise MadGraph5Error("Could not determine the type of integrand to be added for"+
                                     " the contribution definiton:\n%s"%str(contribution_definition.nice_string()))
            return super(ME7Integrand, cls).__new__(target_class, *all_args, **opt)
        else:
            return super(ME7Integrand, cls).__new__(cls, *all_args, **opt)
    
    def __init__(self, model, 
                       run_card,
                       contribution_definition,
                       processes_map,
                       all_MEAccessors,
                       ME7_configuration):
        """ Initializes a generic ME7 integrand defined by the high-level abstract information
        given in arguments, whose explanations are provided in the comments of this initializer's body. """
        
        super(ME7Integrand, self).__init__()
        
        # Keep the initialization inputs to generate the dump
        # Only keep the relevant information, in and set to None the ones that will need to be updated by
        # ME7 interface directly.
        self.initialization_inputs = { 'model'                      : None,
                                       'run_card'                   : None,
                                       'contribution_definition'    : contribution_definition,
                                       'processes_map'              : processes_map,
                                       'all_MEAccessors'            : None,
                                       'ME7_configuration'          : None }
        
        # A ModelReader instantance, initialized with the values of the param_card.dat of this run
        self.model                      = model
        if not isinstance(self.model, model_reader.ModelReader):
            raise MadGraph5Error("The ME7Integrand must be initialized with a ModelReader instance.")
        # A RunCardME7 instance, properly initialized with the values of the run_card.dat of this run
        self.run_card                   = run_card
        # The original ContributionDefinition instance at the origin this integrand 
        self.contribution_definition    = contribution_definition
        # The process map of the Contribution instance at the origin of this integrand.
        # The format is identical to the one generated from the function 'get_process_map' of a contribution.
        self.processes_map              = processes_map
        # An instance of contributions.MEAccessorDict providing access to all ME available as part of this
        # ME7 session.
        self.all_MEAccessors            = all_MEAccessors
        # The option dictionary of ME7
        self.ME7_configuration          = ME7_configuration
        
        all_processes = [p[0] for p in self.processes_map.values()]
        self.masses = self.get_external_masses_for_process(all_processes[0], model=self.model)
        for proc in all_processes[1:]:
            this_proc_masses = self.get_external_masses_for_process(proc, model=self.model)
            if this_proc_masses != self.masses:
                raise MadGraph5Error("A contribution must entail processes with all the same external masses.\n"
                 "This is not the case; process\n%s\nhas masses '%s' while process\n%s\n has masses '%s'."%
                 (all_processes[0].nice_string(), masses, proc.nice_string(), this_proc_masses) )
        self.n_initial = len(self.masses[0])
        self.n_final = len(self.masses[1])
        
        # Always initialize the basic flat PS generator. It can be overwritten later if necessary.
        self.phase_space_generator = phase_space_generators.FlatInvertiblePhasespace(self.masses[0], self.masses[1])
        integrand_dimensions = integrands.DimensionList()
        
        if self.n_initial==1:
            raise InvalidCmd("MadEvent7 does not yet support decay processes.")
        
        # Add the PDF dimensions if necessary
        if self.run_card['lpp1']==self.run_card['lpp2']==1:
            self.collider = 'pp'
            integrand_dimensions.append(integrands.ContinuousDimension('ycms',lower_bound=0.0, upper_bound=1.0))
            # The 2>1 topology requires a special treatment
            if self.n_initial==2 and self.n_final==1:
                self.collider = 'pp_2to1'
            else:
                integrand_dimensions.append(integrands.ContinuousDimension('tau',lower_bound=0.0, upper_bound=1.0)) 
        elif self.run_card['lpp1']==self.run_card['lpp2']==0:
            self.collider = 'll'
        else:
            raise InvalidCmd("MadEvent7 does not support this collider configuration: (lpp1=%d, lpp2=%d)"%
                             (self.run_card['lpp1'], self.run_card['lpp2']))
        
        if self.run_card['ebeam1']!=self.run_card['ebeam2']:
            raise InvalidCmd("For now, MadEvent7 only supports colliders with incoming beams equally energetic.")

        # Add the phase-space dimensions
        integrand_dimensions.extend([ integrands.ContinuousDimension('x_%d'%i,lower_bound=0.0, upper_bound=1.0) 
                                     for i in range(1, self.phase_space_generator.nDimPhaseSpace()+1) ])
        self.set_dimensions(integrand_dimensions)

        self.collider_energy = self.run_card['ebeam1'] + self.run_card['ebeam2']
        
        # Set the seed
        if self.run_card['iseed'] > 0:
            random.seed(self.run_card['iseed'])

        # Initialize the PDF, if necessary
        # Setup the PDF cache
        self.PDF_cache = {}
        self.PDF_cache_entries = []
        if self.run_card['lpp1']==0 and self.run_card['lpp2']==0:
            self.pdf = None 
            self.pdfsets = None
        else:
            if self.run_card['pdlabel'] != 'lhapdf':
                raise InvalidCmd("MadEvent7 does not support built-in PDFs.")
            # Load the PDFs
            if self.ME7_configuration['lhapdf']:
                lhapdf_config = self.ME7_configuration['lhapdf']
            else:
                lhapdf_config = misc.which('lhapdf-config')
            lhapdf = misc.import_python_lhapdf(lhapdf_config)
            if not lhapdf:
                raise MadGraph5Error("The python lhapdf API could not be loaded.")
            # Adjust LHAPDF verbosity to current logger's verbosity
            lhapdf.setVerbosity(1 if logger.level<=logging.DEBUG else 0)
            
            pdfsets_dir = subprocess.Popen([lhapdf_config,'--datadir'],\
                                           stdout=subprocess.PIPE).stdout.read().strip()
            lhapdf_version = subprocess.Popen([lhapdf_config,'--version'],\
                                           stdout=subprocess.PIPE).stdout.read().strip()
            pdf_info = common_run.CommonRunCmd.get_lhapdf_pdfsets_list_static(pdfsets_dir, lhapdf_version)
            lhaid = self.run_card.get_lhapdf_id()
            if lhaid not in pdf_info:
                raise InvalidCmd("Could not find PDF set with lhaid #%d in %s."%(lhaid, pdfsets_dir))  
            self.pdfsets = lhapdf.getPDFSet(pdf_info[lhaid]['filename'])
            # Pick the central PDF for now
            self.pdf = self.pdfsets.mkPDF(0)

        # Initialize access to a mean of running alpha_s
        if self.pdf:
            self.alpha_s_runner = self.pdf.alphasQ
        else:
            # We assume here that the model passed to this integrand for initialization started off
            # with its couplings (i.e. alpha_s) defined for Q**2= MZ**2.
            model_param_dict = model.get('parameter_dict')
            as_running_params = {'aS':0.118, 'mdl_MZ':91.188, 'mdl_MC':1.55, 'mdl_MB':4.7}
            for param in as_running_params:
                if param not in model_param_dict:
                    if param in ['mdl_MC', 'mdl_MB']:
                        # Leave these parameters to their default if not specified in the model
                        continue
                    else:
                        raise InvalidCmd("When not using PDFsets, MadEvent7 requires a model with the"+
                                    " parameter %s to be defined so as to be able to run alpha_S."%param)
                if model_param_dict[param] != 0.:
                    as_running_params[param] = model_param_dict[param]
            # For now always chose to run alpha_S at two loops.
            n_loop_for_as_running = 2
            self.alpha_s_runner = model_reader.Alphas_Runner(as_running_params['aS'], n_loop_for_as_running, 
                      as_running_params['mdl_MZ'], as_running_params['mdl_MC'], as_running_params['mdl_MB'])

    def generate_dump(self):
        """ Generate a serializable dump of self, which can later be used, along with some more 
        information, in initialize_from_dump in order to regenerate the object."""
        
        dump = {'class': self.__class__}
        dump.update(self.initialization_inputs)
        return dump
    
    @classmethod
    def initialize_from_dump(cls, dump, model, run_card, all_MEAccessors, ME7_configuration):
        """ Initialize self from a dump and possibly other information necessary for reconstructing this
        integrand."""
        
        return dump['class'](model, 
                             run_card,                             
                             dump['contribution_definition'],
                             dump['processes_map'],
                             all_MEAccessors,
                             ME7_configuration)
      
    def get_external_masses_for_process(self, process, model):
        """ Returns the tuple:
               ( (initial_mass_value1, ...) , (final_mass_value1, final_mass_value2, final_mass_value3,...)
            for the process in argument
        """

        return ( tuple(model.get_mass(pdg) for pdg in process.get_initial_ids()),
                 tuple(model.get_mass(pdg) for pdg in process.get_final_ids()),
               )
    
    def set_phase_space_generator(self, PS_generator):
        """ Overwrites current phase-space generator."""
        if not isinstance(PS_generator, phase_space_generators.VirtualPhaseSpaceGenerator):
            raise MadGraph5Error("Cannot assign to a MadEvent7 integrand a phase-space generator that "+
                                 " does not inherit from VirtualPhaseSpaceGenerator.")
        if PS_generator.nDimPhaseSpace() != self.phase_space_generator.nDimPhaseSpace():
            raise MadGraph5Error("A MadEvent7 integrand was assigned a phase-space generator with the"+
                                 " wrong number of integration dimensions: %d instead of %d"%
                (PS_generator.nDimPhaseSpace(),self.phase_space_generator.nDimPhaseSpace()))
        self.phase_space_generator = PS_generator
        
    def pass_flavor_blind_cuts(self, PS_point, process_pdgs):
        """ Implementation of a minimal set of isolation cuts. This can be made much nicer in the future and 
        will probably be taken outside of this class so as to ease user-specific cuts, fastjet, etc...
        We consider here a two-level cuts system, this first one of which is flavour blind.
        This is of course not IR safe at this stage!"""

        # These cuts are not allowed to resolve flavour, but only whether a particle is a jet or not
        def is_a_jet(pdg):
            return pdg in range(1,7)+range(-1,-7,-1)+[21]

        logger.debug( "Processing flavor-blind cuts for process %s and PS point:\n%s"%(
                        str(process_pdgs), self.phase_space_generator.nice_momenta_string(
                            PS_point, recompute_mass=True, n_initial=self.phase_space_generator.n_initial) ))

        pt_cut = self.run_card['ptj']
        dr_cut = self.run_card['drjj']
        
        if pt_cut < 0. and dr_cut < 0.:
            return True
        
        # Apply the Ptj cut first
        for i, p in enumerate(PS_point[self.n_initial:]):
            if not is_a_jet(process_pdgs[self.n_initial+i]):
                continue
            logger.debug('p_%i.pt()=%.5e'%((self.n_initial+i),p.pt()))
            if p.pt() < pt_cut:
                return False

        # And then the drjj cut
        for i, p1 in enumerate(PS_point[self.n_initial:]):
            for j, p2 in enumerate(PS_point[self.n_initial+i+1:]):
                if not is_a_jet(process_pdgs[self.n_initial+i]) and\
                   not is_a_jet(process_pdgs[self.n_initial+i+1+j]):
                    continue
                logger.debug('deltaR(p_%i,p_%i)=%.5e'%(
                     self.n_initial+i, self.n_initial+i+1+j, p1.deltaR(p2)))
                if p1.deltaR(p2) < dr_cut:
                    return False

        return True

    def pass_flavor_sensitive_cuts(self, PS_point, flavors):
        """ Implementation of a minimal set of isolation cuts. This can be made much nicer in the future and 
        will probably be taken outside of this class so as to ease user-specific cuts, fastjet, etc...
        We consider here a two-level cuts system, this second one of which is flavour sensitive."""

        logger.debug( "Processing flavor-sensitive cuts for flavors %s and PS point:\n%s"%(
                        str(flavors), self.phase_space_generator.nice_momenta_string(
                            PS_point, recompute_mass=True, n_initial=self.phase_space_generator.n_initial) ))

        # None implemented yet
        return True

    @staticmethod
    def Lambda(s,sqrMA,sqrMB):
        """ Kahlen function."""

        return s**2 + sqrMA**2 + sqrMB**2 - 2.*s*sqrMA - 2.*sqrMB*sqrMA - 2.*s*sqrMB

    def get_scales(self, PS_point):
        """ Returns mu_r, mu_f1, mu_f2 for that PS point."""
        
        if not self.run_card['fixed_ren_scale'] or not self.run_card['fixed_fac_scale']:
            raise InvalidCmd("MadEvent7 only supports fixed mu_f and mu_r scales for now.")
        
        return self.run_card['scale'], self.run_card['dsqrt_q2fact1'], self.run_card['dsqrt_q2fact2'] 

    def get_pdfQ2(self, pdf, pdg, x, scale2):
       
        if pdg not in [21,22] and abs(pdg) not in range(1,7):
            return 1.
                
        if (pdf, pdg, x, scale2) in self.PDF_cache:
            return self.PDF_cache[(pdf, pdg, x, scale2)]
        
        # Call to lhapdf API
        f = pdf.xfxQ2(pdg, x, scale2)/x
        
        # Update the PDF cache
        self.PDF_cache[(pdf, pdg,x,scale2)] = f
        self.PDF_cache_entries.append((pdf, pdg,x,scale2)) 
        if len(self.PDF_cache_entries) > self.PDF_cache_max_size:
            del self.PDF_cache[self.PDF_cache_entries.pop(0)]

        return f 

    def __call__(self, continuous_inputs, discrete_inputs, **opts):
        """ Main function of the integrand, returning the weight to be passed to the integrator."""
        
        # A unique float must be returned
        wgt = 1.0
        
        logger.debug("="*80)       
        logger.debug('Starting a new evaluation of the integrand from contribution:\n%s',
                                                    self.contribution_definition.nice_string())
        
        # Random variables sent
        random_variables    = list(continuous_inputs)
        logger.debug('Random variables received: %s',str(random_variables))        
    
        # Avoid extrema since the phase-space generation algorithm doesn't like it
        random_variables = [min(max(rv,self.epsilon_border),1.-self.epsilon_border) for rv in random_variables]
        
        # Assign variables to their meaning. We could use their name to be more definite,
        # but assuming their position is is fine.
        dimension_names = [dim.name for dim in self.get_dimensions()]
        dimension_name_to_position = dict((name,i) for i, name in enumerate(dimension_names))
        if 'ycms' in dimension_names:
            PDF_ycm = random_variables[dimension_name_to_position['ycms']]
        else:
            PDF_ycm = None
        if 'tau' in dimension_names:
            PDF_tau = random_variables[dimension_name_to_position['tau']]
        else:
            PDF_tau = None
        PS_random_variables  = [rv for i, rv in enumerate(random_variables) if dimension_names[i].startswith('x') ]

        # And the conversion from GeV^-2 to picobarns
        wgt *= 0.389379304e9

        # Now take care of the Phase-space generation:
        
        # Set some defaults for the variables to be set further
        xb_1 = 1.
        xb_2 = 1.
        E_cm = self.collider_energy
        
        # We generate the PDF from two variables \tau = x1*x2 and ycm = 1/2 * log(x1/x2), so that:
        #  x_1 = sqrt(tau) * exp(ycm)
        #  x_2 = sqrt(tau) * exp(-ycm)
        # The jacobian of this transformation is 1.
        
        if self.collider.startswith('pp'):
            
            tot_final_state_masses = sum(self.masses[1])
            if tot_final_state_masses > self.collider_energy:
                raise InvalidCmd("Collider energy is not large enough, there is no phase-space left.")
            
            # Keep a hard cut at 1 GeV, which is the default for absolute_Ecm_min
            tau_min = (max(tot_final_state_masses, self.absolute_Ecm_min)/self.collider_energy)**2
            tau_max = 1.0

            if self.collider == 'pp_2to1':
                # Here tau is fixed by the \delta(xb_1*xb_2*s - m_h**2) which sets tau to 
                PDF_tau = tau_min
                # Account for the \delta(xb_1*xb_2*s - m_h**2) and corresponding y_cm matching to unit volume
                wgt *= (1./self.collider_energy**2)
            else:
                # Rescale tau appropriately
                PDF_tau = tau_min+(tau_max-tau_min)*PDF_tau
                # Including the corresponding Jacobian
                wgt *= (tau_max-tau_min)

            # And we can now rescale ycm appropriately
            ycm_min = 0.5 * math.log(PDF_tau)
            ycm_max = -ycm_min
            PDF_ycm = ycm_min + (ycm_max - ycm_min)*PDF_ycm            
            # and account for the corresponding Jacobina
            wgt *= (ycm_max - ycm_min)

            xb_1 = math.sqrt(PDF_tau) * math.exp(PDF_ycm)
            xb_2 = math.sqrt(PDF_tau) * math.exp(-PDF_ycm)
            E_cm = math.sqrt(xb_1*xb_2)*self.collider_energy

        elif self.collider == 'll':
            xb_1 = 1.
            xb_2 = 1.
            E_cm = self.collider_energy
        else:
            raise MadGraph5Error("MadEvent7 integrand does not yet support collider mode '%s'."%self.collider)

        # Make sure the Bjorken x's are physical:
        if xb_1>1. or xb_2>1.:
            logger.debug('Unphysical configuration: x1, x2 = %.5e, %.5e'%(xb_1, xb_2))
            logger.debug(misc.bcolors.GREEN + 'Returning a weight of 0. for this integrand evaluation.' + misc.bcolors.ENDC)
            return 0.0

        # Include the flux factor
        flux = 1.
        if self.n_initial == 2:
            flux = 1. / (2.*math.sqrt(self.Lambda(E_cm**2, self.masses[0][0], self.masses[0][1])))
        elif self.n_initial == 1:
            flux = 1. / (2.*E_cm)
        flux /= math.pow(2.*math.pi, 3*self.n_final - 4)
        wgt *= flux
        logger.debug("Flux factor: %.5e"%flux)

        # Now generate a PS point
        PS_point, PS_weight = self.phase_space_generator.generateKinematics(E_cm, PS_random_variables)

        # Account for PS weight
        wgt *= PS_weight
        logger.debug("PS_weight: %.5e"%PS_weight)

        ###
        # /!\ WARNING ONLY BOOST TO C.O.M frame at the very end. #
        ###

        # We must now boost the PS point in the lab frame
        if self.n_initial == 2:
            ref_lab = PS_point[0]*xb_1 + PS_point[1]*xb_2
            if ref_lab.rho2() != 0.:
                ref_lab.setMass(ref_lab.calculateMass())
                for p in PS_point:
                    p.boost(ref_lab.boostVector())

        logger.debug("Considering the following PS point:\n%s"%(self.phase_space_generator.nice_momenta_string(
                            PS_point, recompute_mass=True, n_initial=self.phase_space_generator.n_initial) ))

        # Recover scales to be used
        mu_r, mu_f1, mu_f2 = self.get_scales(PS_point)

        # Apply the recomputed alpha_s from the PDFs or a dedicated model_reader.Alphas_Runner
        alpha_s = self.alpha_s_runner(mu_r**2)
        # Notice here that we do *not* update all the couplings/parameters dependent on this new value of mu_r / alpha_s
        # We reset only mu_r and alpha_s since this is the only thing our integrands directly depend on so far (along with
        # the masses and widths which are of course alpha_s independent.)
        # Indeed the matrix element take their couplings and others directly from the fortran exported version (i.e via f2py).
        # Also it can be quite time consuming to update the whole set of dependent parameters that are python expression, so
        # it is best to avoid it if possible.
        model_param_dict = self.model.get('parameter_dict')
        model_param_dict['aS'] = alpha_s
        if 'MU_R' in model_param_dict:
            model_param_dict['MU_R'] = mu_r
            
        # Now loop over processes
        total_wgt = 0.
        for process_hash, (process, mapped_processes) in self.processes_map.items():
            logger.debug('Now considering the process group from %s.'%process.nice_string())
            
            this_process_wgt = wgt
            process_pdgs = tuple(process.get_initial_ids()+process.get_final_ids())

            # Apply flavor blind cuts
            if not self.pass_flavor_blind_cuts(PS_point, process_pdgs):
                logger.debug('Event failed the flavour_blind generation-level cuts.')
                logger.debug(misc.bcolors.GREEN + 'Returning a weight of 0. for this integrand evaluation.' + misc.bcolors.ENDC)            
                return 0.0
            
            all_process_pdgs = [process_pdgs,] + \
                               [tuple(proc.get_initial_ids()+proc.get_final_ids()) for proc in mapped_processes]

            if self.collider.startswith('pp'):
                logger.debug("Bjorken x's x1, x2, sqrt(x1*x2*s): %.5e, %.5e. %.5e"%(
                                                            xb_1, xb_2, math.sqrt(self.collider_energy**2*xb_1*xb_2)))
            proc_PDFs_weights = []
            for proc_pdgs in all_process_pdgs:            
                if self.collider.startswith('pp'):            
                    PDF1 = self.get_pdfQ2(self.pdf, proc_pdgs[0], xb_1, mu_f1**2)
                    PDF2 = self.get_pdfQ2(self.pdf, proc_pdgs[1], xb_2, mu_f2**2)
                    proc_PDFs_weights.append(PDF1*PDF2)
                    logger.debug("PDF(x1, %d), PDF(x2, %d) = %.5e, %.5e"%(proc_pdgs[0], proc_pdgs[1], PDF1, PDF2))
                else:
                    proc_PDFs_weights.append(1.)

            # Pick a flavor combination
            abs_pdf_wgts_sum = sum(abs(pdf_wgt) for pdf_wgt in proc_PDFs_weights)
            rv_flavor = random.random()*abs_pdf_wgts_sum
            index_selected=0
            running_sum = abs(proc_PDFs_weights[index_selected])
            while rv_flavor > running_sum:
                index_selected +=1
                running_sum += abs(proc_PDFs_weights[index_selected])

            selected_flavors = all_process_pdgs[index_selected]
            logger.debug('Selected flavor combination : %s'%str(selected_flavors))

            # Apply flavor sensitive cuts
            if not self.pass_flavor_sensitive_cuts(PS_point, selected_flavors):
                logger.debug('Event failed the flavour_sensitive generation-level cuts.')
                logger.debug(misc.bcolors.GREEN + 'Returning a weight of 0. for this integrand evaluation.' + misc.bcolors.ENDC)            
                return 0.0      
            
            # Apply PDF weight
            this_flavor_wgt = this_process_wgt*proc_PDFs_weights[index_selected]
            this_process_wgt *= sum(proc_PDFs_weights)
    
            # Finally include the short-distance weight
            selected_flavors = ( tuple(selected_flavors[i] for i in range(self.n_initial)),
                                 tuple(selected_flavors[self.n_initial+i] for i in range(self.n_final)))

            sigma_wgt = self.sigma(PS_point, process, selected_flavors, this_flavor_wgt, mu_r, mu_f1, mu_f2)
            logger.debug('Short-distance sigma weight for this subprocess: %.5e'%sigma_wgt)        
            this_process_wgt *= sigma_wgt
            
            # Accumulate this process weight
            total_wgt += this_process_wgt

        # Now finally return the total weight for this contribution
        logger.debug(misc.bcolors.GREEN + "Final weight returned: %.5e"%total_wgt + misc.bcolors.ENDC)
        logger.debug("="*80)
        return total_wgt
    
    def sigma(self, PS_point, process, flavors, flavor_wgt, mu_r, mu_f1, mu_f2):
        """ 
        This is the core function of the integrand where the short-distance objects like the matrix elements,
        the counterterms, the mappings, etc.. will be evaluated.
        The process has a lot of extra meta information that sigma may consider, but flavors provides what flavour 
        assignment should be used for some of the flavour dependent steps (like calling the cuts for instance).
        Notice that this function may be called several times with the exact same arguments but different 
        processes / flavour, so some caching may be in order. 
        The output of sigma, sigma_wgt, should be the weight summed over the various contributions building sigma
        (ME's and counterterms).
        Finally, this function is also responsible for calling the observable, using flavor_wgt*sigma_wgt.
        It can potentially call it several times, for the various counterterms.
        """

        sigma_wgt = 1.

        alpha_s = self.model.get('parameter_dict')['aS']

        ME_evaluation = self.all_MEAccessors(process, PS_point, alpha_s, mu_r, pdgs=flavors)

        # Notice here that the most general call to the Matrix Element would be:
        #
        # ME_evaluation = self.all_MEAccessors(process, PS_point, alpha_s, mu_r, pdgs=flavors,
        #       squared_orders = {...}, 
        #       spin_correlation = [...],
        #       color_connection = [...], 
        #       hel_config = [...] )
        #
        # One can read the details of the format for each of these options in
        # the documentation of the function contributions.VirtualMEAccessor.apply_permutations
        
        # Also, here is the more pedantic way of obtaining an ME evaluation:
        ## process_key = contributions.ProcessKey(process=process, pdgs=flavors)
        ## ME_accessor, call_key = self.all_MEAccessors.get_MEAccessor(process_key, pdgs=flavors)
        ## call_key['squared_orders'] = {...}
        ## call_key['spin_correlation'] = [...]
        ## call_key['color_connection'] = [...]
        ## call_key['color_connection'] = [...]
        ## ...
        ## ME_evaluation = ME_accessor(PS_point, alpha_s, mu_r, **call_key)

        sigma_wgt *= ME_evaluation['finite']
        
        if self.apply_observables:
            data_for_observables = {'PS_point': PS_point, 'flavors' : flavors}
            self.observable_list.apply_observables(sigma_wgt*flavor_wgt, data_for_observables)

        return sigma_wgt

# Daughter classes for the various type of Integrand/contributions which will surely require a redefinition of sigma()
# and possibly other functions.
class ME7Integrand_B(ME7Integrand):
    """ ME7Integrand for the computation of a Born type of contribution."""
    def sigma(self, PS_point, process, flavors, flavor_wgt, mu_r, mu_f1, mu_f2, *args, **opts):
        return super(ME7Integrand_B, self).sigma(PS_point, process, flavors, flavor_wgt, mu_r, mu_f1, mu_f2, *args, **opts)

class ME7Integrand_LI(ME7Integrand):
    """ ME7Integrand for the computation of a Loop-Induced Born type of contribution."""
    def sigma(self, PS_point, process, flavors, flavor_wgt, mu_r, mu_f1, mu_f2, *args, **opts):
        return super(ME7Integrand_LI, self).sigma(PS_point, process, flavors, flavor_wgt, mu_r, mu_f1, mu_f2, *args, **opts)
        
class ME7Integrand_V(ME7Integrand):
    """ ME7Integrand for the computation of a one-loop virtual type of contribution."""
    def sigma(self, PS_point, process, flavors, flavor_wgt, mu_r, mu_f1, mu_f2, *args, **opts):
        return super(ME7Integrand_V, self).sigma(PS_point, process, flavors, flavor_wgt, mu_r, mu_f1, mu_f2, *args, **opts)

class ME7Integrand_R(ME7Integrand):
    """ ME7Integrand for the computation of a single real-emission type of contribution."""
    def sigma(self, PS_point, process, flavors, flavor_wgt, mu_r, mu_f1, mu_f2, *args, **opts):
        return super(ME7Integrand_R, self).sigma(PS_point, process, flavors, flavor_wgt, mu_r, mu_f1, mu_f2, *args, **opts)

class ME7Integrand_RR(ME7Integrand):
    """ ME7Integrand for the computation of a double real-emission type of contribution."""
    def sigma(self, PS_point, process, flavors, flavor_wgt, mu_r, mu_f1, mu_f2, *args, **opts):
        return super(ME7Integrand_RR, self).sigma(PS_point, process, flavors, flavor_wgt, mu_r, mu_f1, mu_f2, *args, **opts)
