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
logger = logging.getLogger('madevent7') # -> stdout
logger_stderr = logging.getLogger('madevent7.stderr') # ->stderr

try:
    import madgraph
except ImportError: 
    # import from madevent directory
    MADEVENT7 = True
    import internal.extended_cmd as cmd
    import internal.common_run_interface as common_run
    import internal.madevent_interface as madevent_interface
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
    import madgraph.interface.madevent_interface as madevent_interface
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
    
    def parse_test_IR_limits(self, args):
        """ Parses the options specified to the test_limit command."""
        
        # None means unspecified, therefore considering all types.
        testlimits_options = {'correction_order'        : None,
                              'limit_type'              : None,
                              'process'                 : {'in_pdgs'  : None,
                                                           'out_pdgs' : None,
                                                           'n_loops'  : None},
                              'seed'                    : None,
                              'n_steps'                 : 10,
                              'min_scaling_variable'    : 1.0e-6,
                              'acceptance_threshold'    : 1.0e-6
                             }
 
        # Group arguments in between the '--' specifiers.
        # For example, it will group '--process=p' 'p' '>' 'd' 'd~' 'z' into '--process=p p > d d~ z'.
        opt_args = []
        new_args = []
        for arg in args:
            if arg.startswith('--'):
                opt_args.append(arg)
            elif len(opt_args)>0:
                opt_args[-1] += ' %s'%arg
            else:
                new_args.append(arg)
        
        for arg in opt_args:
            try:
                key, value = arg.split('=',1)
                value = value.strip()
            except ValueError:
                key = arg
                value = None
            
            if key == '--correction_order':
                if value.upper() not in ['NLO','NNLO','NNNLO']:
                    raise InvalidCmd("'%s' is not a valid option for '%s'"%(value, key))
                testlimits_options['correction_order'] = value.upper()
            elif key in ['--n_steps']:
                try:
                    testlimits_options[key[2:]] = int(value)
                    if int(value)<2:
                        raise ValueError  
                except ValueError:
                    raise InvalidCmd("'%s' is not a valid integer for option '%s'"%(value, key))
            elif key == '--n_loops':
                try:
                    testlimits_options['process'][key[2:]] = int(value)
                    if int(value)<0:
                        raise ValueError
                except ValueError:
                    raise InvalidCmd("'%s' is not a valid integer for option '%s'"%(value, key))
            elif key in ['--min_scaling_variable', '--acceptance_threshold']:
                try:
                    testlimits_options[key[2:]] = float(value)                  
                except ValueError:
                    raise InvalidCmd("'%s' is not a valid float for option '%s'"%(value, key))                  
            elif key in '--limit_type':
                if not isinstance(value, str):
                    raise InvalidCmd("'%s' is not a valid option for '%s'"%(value, key))
                if value.lower() == 'soft':
                    testlimits_options['limit_type'] = re.compile(r'.*S.*')                    
                elif value.lower() == 'collinear':
                    testlimits_options['limit_type'] = re.compile(r'.*C.*')
                elif value.lower() == 'all':
                    testlimits_options['limit_type'] = re.compile(r'.*')
                else:
                    if any(value.startswith(start) for start in ['r"',"r'"]):
                        testlimits_options['limit_type'] = re.compile(value)
                    else:
                        # If the specified re was not explicitly made a raw string, then we take the 
                        # liberty here of escaping the parenthesis since this is presumably what the
                        # user expects.
                        testlimits_options['limit_type'] = re.compile(value.replace('(','\(').replace(')','\)'))
                        
            elif key == '--seed':
                try:
                    testlimits_options['seed'] = int(value)
                except ValueError:
                    raise InvalidCmd("Cannot set '%s' option to '%s'."%(key, value))
            elif key=='--process':
                def get_particle_pdg(name):
                    part = self.model.get_particle(name)
                    part.set('is_part', name!=part.get('antiname'))
                    return part.get_pdg_code()
                initial, final = value.split('>',1)
                initial = initial.strip()
                final = final.strip()
                initial_pdgs = []
                for parts in initial.split(' '):
                    multi_part = []
                    for part in parts.split('|'):
                        try:
                            multi_part.append(get_particle_pdg(part))
                        except:
                            raise InvalidCmd("Particle '%s' not recognized in current model."%part)
                    initial_pdgs.append(tuple(multi_part))

                final_pdgs = []
                for parts in final.split(' '):
                    multi_part = []
                    for part in parts.split('|'):
                        try:
                            multi_part.append(get_particle_pdg(part))
                        except:
                            raise IvalidCmd("Particle '%s' not recognized in current model."%part)
                    final_pdgs.append(tuple(multi_part))
                    
                testlimits_options['process']['in_pdgs'] = tuple(initial_pdgs)
                testlimits_options['process']['out_pdgs'] = tuple(final_pdgs)
            else:
                raise InvalidCmd("Option '%s' for the test_limits command not recognized."%key)        
        
        return new_args, testlimits_options
    
    def parse_launch(self, args):
        """ Parses the argument of the launch command."""

        launch_options = {'integrator': 'VEGAS3',
                          'n_points': None,
                          'n_iterations':None,
                          'verbosity':1,
                          'refresh_filters':'auto',
                          'compile':'auto'}        
        
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
            elif key in ['--refresh_filters','--compile']:
                available_modes = ['auto','never','always']
                if value is None:
                    mode = 'always'
                else:
                    mode = str(value).lower()
                if mode not in available_modes:
                    raise InvalidCmd("Value '%s' not valid for option '%s'."%(value,key))
                launch_options[key[2:]] = mode
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
        self.prompt = "ME7::%s > "%os.path.basename(pjoin(self.me_dir))
        
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

        self.complex_mass_scheme = ME7_dump['model_with_CMS']
                
        # We might want to recover whether prefix was used when importing the model and whether
        # the MG5 name conventions was used. But this is a detail that can easily be fixed later.
        self.model = import_ufo.import_model(
            pjoin(self.me_dir,'Source',ME7_dump['model_name']), prefix=True,
            complex_mass_scheme = ME7_dump['model_with_CMS'] )
        self.model.pass_particles_name_in_mg_default()
        
        self.model.set_parameters_and_couplings(
                param_card = pjoin(self.me_dir,'Cards','param_card.dat'), 
                scale=self.run_card['scale'], 
                complex_mass_scheme=self.complex_mass_scheme)

        self.all_MEAccessors = ME7_dump['all_MEAccessors']['class'].initialize_from_dump(
                                ME7_dump['all_MEAccessors'], root_path = self.me_dir, model=self.model)
        self.all_integrands = [integrand_dump['class'].initialize_from_dump(
            integrand_dump, self.model, self.run_card, self.all_MEAccessors, self.options
                                                    ) for integrand_dump in ME7_dump['all_integrands']]
        
        self.mode = self.get_maximum_overall_correction_order()
        self.prompt = "ME7@%s::%s > "%(self.mode, os.path.basename(pjoin(self.me_dir)))

        self.n_initial = ME7_dump['n_initial']
        
    def get_maximum_overall_correction_order(self):
        """ From investigating the integrands, this function derives what is the maximum correction
        order included in this ME7 session."""
        max_order = 'LO'
        for integrand in self.all_integrands:
            max_integrand_order = integrand.contribution_definition.overall_correction_order 
            if max_integrand_order.count('N') > max_order.count('N'):
                max_order = max_integrand_order
        return max_integrand_order
        
    def synchronize(self, **opts):
        """ Re-compile all necessary resources and sync integrands with the cards and model"""
        
        logger.info("Synchronizing MadEvent7 internal status with cards and matrix elements source codes...")
    
        self.run_card = banner_mod.RunCardME7(pjoin(self.me_dir,'Cards','run_card.dat'))
        self.model.set_parameters_and_couplings(
            param_card = pjoin(self.me_dir,'Cards','param_card.dat'), 
            scale=self.run_card['scale'], 
            complex_mass_scheme=self.complex_mass_scheme)
        
        for integrand in self.all_integrands:
            integrand.synchronize(self.model, self.run_card, self.options)
        
        # Try and import some options from those provided to this function
        sync_options = {'refresh_filters':'auto', 'compile':'auto', 'model':self.model}
        for key in sync_options:
            try:
                sync_options[key] = opts[key]
            except KeyError:
                pass
        self.all_MEAccessors.synchronize(ME7_options=self.options, **sync_options)
        
    def do_launch(self, line, *args, **opt):
        """Main command, starts the cross-section computation. Very basic setup for now.
        We will eventually want to have all of these meta-data controllable via user commands
        in the interface.
        This is super naive and only for illustrative purposes for now."""
        
        args = self.split_arg(line)
        
        launch_options = self.parse_launch(args)

        # In principle we want to start by recompiling the process output so as to make sure
        # that everything is up to date.
        self.synchronize(**launch_options)

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
            
        logger.info("="*100)        
        logger.info('{:^100}'.format("Starting integration, lay down and enjoy..."),'$MG:color:GREEN')
        logger.info("="*100)

        with misc.MuteLogger(['contributions','madevent7','madevent7.stderr'],[logger_level,logger_level,logger_level]):
            xsec, error = self.integrator.integrate()
            
        logger.info("="*100)
        logger.info('{:^100}'.format("Cross-section with integrator '%s':"%self.integrator.get_name()),'$MG:color:GREEN')
        logger.info('{:^100}'.format("%.5e +/- %.2e [pb]"%(xsec, error)),'$MG:color:BLUE')
        logger.info("="*100+"\n")
    
    def do_test_IR_limits(self, line, *args, **opt):
        """This function test that local subtraction counterterms match with the actual matrix element in the IR limit."""
    
        args = self.split_arg(line)
        args, testlimits_options = self.parse_test_IR_limits(args)
        
        if testlimits_options['correction_order'] is None:
            # If not defined, automatically assign correction_order to the highest correction considered.
            testlimits_options['correction_order'] = self.mode

        for integrand in self.all_integrands:
            if not hasattr(integrand, 'test_IR_limits'):
                continue
            logger.debug('Now testing IR limits of the following integrand:\n%s'%(integrand.nice_string()))
            integrand.test_IR_limits(test_options = testlimits_options)

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
    pass
    
#===============================================================================
# ME7Integrand
#===============================================================================

class ME7Integrand(integrands.VirtualIntegrand):
    """ Specialization for multi-purpose integration with ME7."""
    
    # Maximum size of the cache for PDF calls
    PDF_cache_max_size = 1000
    
    def __new__(cls, model, 
                     run_card,
                     contribution_definition,
                     processes_map,
                     all_MEAccessors,
                     ME7_configuration, **opt):
        all_args = [model, run_card, contribution_definition, processes_map,
                    all_MEAccessors, ME7_configuration ]
        if cls is ME7Integrand:
            target_type = 'Unknown'
            if contribution_definition.correction_order == 'LO':
                if contribution_definition.n_loops == 0 and \
                  contribution_definition.n_unresolved_particles == 0:
                   target_type = 'Born'
                elif contribution_definition.n_loops == 1 and \
                  contribution_definition.n_unresolved_particles == 0:
                   target_type = 'LoopInduced_Born'
            elif contribution_definition.correction_order == 'NLO':
                if contribution_definition.n_loops == 1 and \
                   contribution_definition.n_unresolved_particles == 0:
                    target_type = 'Virtual'
                elif contribution_definition.n_loops == 0 and \
                     contribution_definition.n_unresolved_particles == 1:
                    target_type = 'SingleReals'
            elif contribution_definition.correction_order == 'NNLO':
                if contribution_definition.n_loops == 0 and \
                   contribution_definition.n_unresolved_particles == 2:
                    target_type = 'DoubleReals'            
                else:
                    raise MadGraph5Error("Some NNLO type of integrands are not implemented yet.")    
            else:
                target_type = 'Unknown'
            target_class = ME7Integrand_classes_map[target_type]
            if not target_class:
                raise MadGraph5Error("Could not determine the class of integrand of type '%s' to be added for"%target_type+
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
                                       'ME7_configuration'          : None,
                                       'options'                    : {} }
        
        # The original ContributionDefinition instance at the origin this integrand 
        self.contribution_definition    = contribution_definition
        # The process map of the Contribution instance at the origin of this integrand.
        # The format is identical to the one generated from the function 'get_process_map' of a contribution.
        self.processes_map              = processes_map
        # An instance of contributions.MEAccessorDict providing access to all ME available as part of this
        # ME7 session.
        self.all_MEAccessors            = all_MEAccessors
        
        # Update and define many properties of self based on the provided run-card and model.
        self.synchronize(model, run_card, ME7_configuration)

    def nice_string(self):
        """ For now simply use the contribution_definition and class name for a nice readable representation."""
        
        res = []
        res.append("Instance of class '%s', with the following contribution definition:"%(self.__class__.__name__))
        res.append('\n'.join(' > %s'%line for line in self.contribution_definition.nice_string().split('\n')))
        return '\n'.join(res)

    def synchronize(self, model, run_card, ME7_configuration):
        """ Synchronize this integrand with the most recent run_card and model."""

        # The option dictionary of ME7
        self.ME7_configuration          = ME7_configuration
        
        # A ModelReader instance, initialized with the values of the param_card.dat of this run
        self.model                      = model
        if not isinstance(self.model, model_reader.ModelReader):
            raise MadGraph5Error("The ME7Integrand must be initialized with a ModelReader instance.")

        # A RunCardME7 instance, properly initialized with the values of the run_card.dat of this run
        self.run_card                   = run_card
        
        # Set external masses
        all_processes = [p[0] for p in self.processes_map.values()]
        self.masses = all_processes[0].get_external_masses(self.model)
        for proc in all_processes[1:]:
            this_proc_masses = proc.get_external_masses(self.model)
            if this_proc_masses != self.masses:
                raise MadGraph5Error("A contribution must entail processes with all the same external masses.\n"
                 "This is not the case; process\n%s\nhas masses '%s' while process\n%s\n has masses '%s'."%
                 (all_processes[0].nice_string(), self.masses, proc.nice_string(), this_proc_masses) )
        self.n_initial = len(self.masses[0])
        self.n_final = len(self.masses[1])
        
        if self.n_initial==1:
            raise InvalidCmd("MadEvent7 does not yet support decay processes.")
        
        if not (self.run_card['lpp1']==self.run_card['lpp2']==1) and \
           not (self.run_card['lpp1']==self.run_card['lpp2']==0):
            raise InvalidCmd("MadEvent7 does not support the following collider mode yet (%d,%d)."%\
                                                            (self.run_card['lpp1'], self.run_card['lpp2']))
        
        # Always initialize the basic flat PS generator. It can be overwritten later if necessary.
        self.phase_space_generator = phase_space_generators.FlatInvertiblePhasespace(
            self.masses[0], self.masses[1], 
            beam_Es    = (self.run_card['ebeam1'], self.run_card['ebeam2']),
            beam_types = (self.run_card['lpp1'], self.run_card['lpp2']),
        )

        # Add a copy of the PS generator dimensions here.
        # Notice however that we could add more dimensions pertaining to this integrand only, and PS generation.
        # This is in particular true for discrete integration dimension like sectors, helicities, etc... 
        self.set_dimensions(integrands.DimensionList(self.phase_space_generator.dimensions))
        self.dim_ordered_names = [d.name for d in self.get_dimensions()]
        self.dim_name_to_position = dict((name,i) for i, name in enumerate(self.dim_ordered_names))
        self.position_to_dim_name = dict((v,k) for (k,v) in self.dim_name_to_position.items())

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
            lhapdf.pathsPrepend(pdfsets_dir)
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
                             ME7_configuration,
                             **dump['options'])
    
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

        if __debug__: logger.debug( "Processing flavor-blind cuts for process %s and PS point:\n%s"%(
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
            if __debug__: logger.debug('p_%i.pt()=%.5e'%((self.n_initial+i),p.pt()))
            if p.pt() < pt_cut:
                return False

        # And then the drjj cut
        for i, p1 in enumerate(PS_point[self.n_initial:]):
            for j, p2 in enumerate(PS_point[self.n_initial+i+1:]):
                if not is_a_jet(process_pdgs[self.n_initial+i]) and\
                   not is_a_jet(process_pdgs[self.n_initial+i+1+j]):
                    continue
                if __debug__: logger.debug('deltaR(p_%i,p_%i)=%.5e'%(
                     self.n_initial+i, self.n_initial+i+1+j, p1.deltaR(p2)))
                if p1.deltaR(p2) < dr_cut:
                    return False

        return True

    def pass_flavor_sensitive_cuts(self, PS_point, flavors):
        """ Implementation of a minimal set of isolation cuts. This can be made much nicer in the future and 
        will probably be taken outside of this class so as to ease user-specific cuts, fastjet, etc...
        We consider here a two-level cuts system, this second one of which is flavour sensitive."""

        if __debug__: logger.debug( "Processing flavor-sensitive cuts for flavors %s and PS point:\n%s"%(
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
        # And the conversion from GeV^-2 to picobarns
        wgt *= 0.389379304e9
        
        if __debug__: logger.debug("="*80)       
        if __debug__: logger.debug('Starting a new evaluation of the integrand from contribution:\n%s',
                                                    self.contribution_definition.nice_string())
        
        # Random variables sent
        random_variables    = list(continuous_inputs)
        if __debug__: logger.debug('Random variables received: %s',str(random_variables))        
    
        # Now assign the variables pertaining to PS generations
        PS_random_variables = [random_variables[self.dim_name_to_position[name]] for name in 
                                                self.phase_space_generator.dim_ordered_names]
        
        PS_point, PS_weight, xb_1, xb_2 = self.phase_space_generator.get_PS_point(PS_random_variables)
        E_cm = math.sqrt(xb_1*xb_2)*self.collider_energy
        
        if PS_point is None:
            if __debug__:
                if xb_1 > 1. or xb_2 > 1.:
                    logger.debug('Unphysical configuration: x1, x2 = %.5e, %.5e'%(xb_1, xb_2))
                    logger.debug(misc.bcolors.GREEN + 'Returning a weight of 0. for this integrand evaluation.' + misc.bcolors.ENDC)
                else:
                    logger.debug('Phase-space generation failed.')
                    logger.debug(misc.bcolors.GREEN + 'Returning a weight of 0. for this integrand evaluation.' + misc.bcolors.ENDC)
            return 0.0
    
        ###
        # /!\ WARNING One typically only boosts to the C.O.M frame at the very end. But we do it here nonetheless. Can easily be changed.
        ###
        self.phase_space_generator.boost_to_COM_frame(PS_point, xb_1, xb_2)
                
        if __debug__: logger.debug("Considering the following PS point:\n%s"%(self.phase_space_generator.nice_momenta_string(
                            PS_point, recompute_mass=True, n_initial=self.phase_space_generator.n_initial) ))
        
        # Account for PS weight
        wgt *= PS_weight
        if __debug__: logger.debug("PS_weight: %.5e"%PS_weight)
        
        # Include the flux factor
        flux = 1.
        if self.n_initial == 2:
            flux = 1. / (2.*math.sqrt(self.Lambda(E_cm**2, self.masses[0][0], self.masses[0][1])))
        elif self.n_initial == 1:
            flux = 1. / (2.*E_cm)
        flux /= math.pow(2.*math.pi, 3*self.n_final - 4)
        wgt *= flux
        if __debug__: logger.debug("Flux factor: %.5e"%flux)

        # Recover scales to be used
        mu_r, mu_f1, mu_f2 = self.get_scales(PS_point)

        # Apply the recomputed alpha_s from the PDFs or a dedicated model_reader.Alphas_Runner
        alpha_s = self.alpha_s_runner(mu_r)
        
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
            if __debug__: logger.debug('Now considering the process group from %s.'%process.nice_string())
            
            this_process_wgt = wgt
            process_pdgs = tuple(process.get_initial_ids()+process.get_final_ids())

            # Apply flavor blind cuts
            if not self.pass_flavor_blind_cuts(PS_point, process_pdgs):
                if __debug__: logger.debug('Event failed the flavour_blind generation-level cuts.')
                if __debug__: logger.debug(misc.bcolors.GREEN + 'Returning a weight of 0. for this integrand evaluation.' + misc.bcolors.ENDC)            
                return 0.0
            
            all_processes = [process,]+mapped_processes
            all_process_pdgs = [ tuple(proc.get_initial_ids()+proc.get_final_ids()) for proc in all_processes]
            # Add mirror processes if present
            all_process_pdgs.extend([ tuple(proc.get_initial_ids()[::-1]+proc.get_final_ids()) for proc in all_processes 
                                                                                        if proc.get('has_mirror_process')])

            if abs(self.run_card['lpp1'])==abs(self.run_card['lpp2'])==1:
                if __debug__: logger.debug("Bjorken x's x1, x2, sqrt(x1*x2*s): %.5e, %.5e. %.5e"%(
                                                            xb_1, xb_2, math.sqrt(self.collider_energy**2*xb_1*xb_2)))
            proc_PDFs_weights = []
            for proc_pdgs in all_process_pdgs:            
                if self.run_card['lpp1']==self.run_card['lpp2']==1:            
                    PDF1 = self.get_pdfQ2(self.pdf, proc_pdgs[0], xb_1, mu_f1**2)
                    PDF2 = self.get_pdfQ2(self.pdf, proc_pdgs[1], xb_2, mu_f2**2)
                    proc_PDFs_weights.append(PDF1*PDF2)
                    if __debug__: logger.debug("PDF(x1, %d), PDF(x2, %d) = %.5e, %.5e"%(proc_pdgs[0], proc_pdgs[1], PDF1, PDF2))
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

            if __debug__: logger.debug('Selected flavor combination : %s'%str(selected_flavors))

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
            if __debug__: logger.debug('Short-distance sigma weight for this subprocess: %.5e'%sigma_wgt)        
            this_process_wgt *= sigma_wgt
            
            # Accumulate this process weight
            total_wgt += this_process_wgt

        # Now finally return the total weight for this contribution
        if __debug__: logger.debug(misc.bcolors.GREEN + "Final weight returned: %.5e"%total_wgt + misc.bcolors.ENDC)
        if __debug__: logger.debug("="*80)
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

        # Access to the matrix element. ME_result is an instance of a subclassed dictionary which includes
        # all independent results available as of now for the particular arguments specified.

#        start = time.time() #TOBECOMMENTED        
        ME_evaluation, all_results = self.all_MEAccessors(process, PS_point, alpha_s, mu_r, pdgs=flavors)
#        misc.sprint("The complete ME call took %10g"%(time.time()-start)) #TOBECOMMENTED
        
        ##
        ## Notice here that the most general call to the Matrix Element would be:
        ##
        ## Example below applies for generate g d > d d d~ QCD^2<=99 QED^2<=99 --LO
        ##
        
        #hel_config = tuple((-1 if fl<0 else 1) for fl in list(flavors[0])+list(flavors[1]))
        #ME_evaluation, ME_result = self.all_MEAccessors(
        #       process, PS_point, alpha_s, mu_r, pdgs=flavors,
        #       squared_orders    = {'QED':2,'QCD':4},
        #       color_correlation = [( 1, 2 )],
        #       spin_correlation  = [( 1, ((1.0,2.0,3.0,4.0), (5.0,6.0,7.0,8.0), (9.0,10.0,11.0,12.0)) )], 
        #       hel_config        = hel_config
        #)
        
        ## Nicely printout the results generated
        #misc.sprint(process.nice_string())
        #misc.sprint(' Flavors: ',flavors)
        #misc.sprint('PS point:\n', self.phase_space_generator.nice_momenta_string(
        #            PS_point, recompute_mass=True, n_initial=self.phase_space_generator.n_initial))
        #misc.sprint('Results pertaining to the specified options:\n'+str(ME_evaluation))
        #misc.sprint('All results generated along with this ME call:\n'+str(all_results))
        ## One can read the details of the format for each of these options in
        ## the documentation of the function contributions.VirtualMEAccessor.apply_permutations
        ## Also, here is the more pedantic way of obtaining an ME evaluation:
        # process_key = contributions.ProcessKey(process=process, pdgs=flavors)
        # ME_accessor, call_key = self.all_MEAccessors.get_MEAccessor(process_key, pdgs=flavors)
        # call_key['squared_orders'] = {...}
        # call_key['spin_correlation'] = [...]
        # call_key['color_correlation'] = [...]
        # ...
        # ME_evaluation = ME_accessor(PS_point, alpha_s, mu_r, **call_key)
        
        ## To debug it is useful to hard-stop the code in a unique noticeable way with a syntax error.
        #stop

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

class ME7Integrand_LIB(ME7Integrand):
    """ ME7Integrand for the computation of a Loop-Induced Born type of contribution."""
    def sigma(self, PS_point, process, flavors, flavor_wgt, mu_r, mu_f1, mu_f2, *args, **opts):
        return super(ME7Integrand_LIB, self).sigma(PS_point, process, flavors, flavor_wgt, mu_r, mu_f1, mu_f2, *args, **opts)
        
class ME7Integrand_V(ME7Integrand):
    """ ME7Integrand for the computation of a one-loop virtual type of contribution."""
    def sigma(self, PS_point, process, flavors, flavor_wgt, mu_r, mu_f1, mu_f2, *args, **opts):
        """ Overloading of the sigma function from ME7Integrand to include necessary additional contributions. """
        
        ret_value = super(ME7Integrand_V, self).sigma(PS_point, process, flavors, flavor_wgt, mu_r, mu_f1, mu_f2, *args, **opts)
        
        return ret_value
        ##
        ## This is just an example to test access to the virtuals
        ##
        ## The example below will work with the process 
        ##     generate u u~ > d d~ QCD^2<=99 QED^2<=99 --NLO=QCD
        ## with NLO correlators forced-in event for the loop.
        
        alpha_s = self.model.get('parameter_dict')['aS']
        hel_config = tuple((-1 if fl<0 else 1) for fl in list(flavors[0])+list(flavors[1]))
        ME_evaluation, all_results = self.all_MEAccessors(
               process, PS_point, alpha_s, mu_r, pdgs=flavors,
               squared_orders    = {'QCD':4,'QED':2},
               color_correlation = [(1, 2), (1,4)],
               spin_correlation  = [( 1, ((1.0,2.0,3.0,4.0), (5.0,6.0,7.0,8.0), (9.0,10.0,11.0,12.0)) )], 
               hel_config        = hel_config 
        )

        ## Nicely print out the results generated.
        misc.sprint(process.nice_string())
        misc.sprint(' Flavors: ',flavors)
        misc.sprint('PS point:\n', self.phase_space_generator.nice_momenta_string(
                    PS_point, recompute_mass=True, n_initial=self.phase_space_generator.n_initial))
        misc.sprint('Results pertaining to the specified options:\n'+str(ME_evaluation))
        misc.sprint('All results generated along with this ME call:\n'+str(all_results))
        
        ## To debug it is useful to hard-stop the code in a unique noticeable way with a syntax error.
        stop
        
        return ret_value

class ME7Integrand_R(ME7Integrand):
    """ ME7Integrand for the computation of a single real-emission type of contribution."""
    
    MappingWalkerType = 'CataniSeymour'
    
    def __init__(self, *args, **opts):
        """ Initialize a real-emission type of integrand, adding additional relevant attributes."""
        
        try:  
            self.counterterms = opts.pop('counterterms')
        except KeyError:
            raise MadEvent7Error("Constructor of class ME7Integrand_R requires the option "+
                                                            "'counterterms' to be specified.")

        super(ME7Integrand_R, self).__init__(*args, **opts)
        self.initialization_inputs['options']['counterterms'] = self.counterterms
    
        # Initialize a general mapping walker capapble of handling all relevant limits for this integrand.
        self.mapper = phase_space_generators.VirtualWalker(map_type=self.MappingWalkerType, model=self.model)
        
    def evaluate_counterterm(self, counterterm, PS_point, hel_config=None):
        """ Evaluates the specified counterterm for the specified PS point."""
        
        # Retrieve some possibly relevant model parameters
        alpha_s = self.model.get('parameter_dict')['aS']
        mu_r = self.model.get('parameter_dict')['MU_R']    
        
        # Now call the mapper to walk through the counterterm structure and return the list of currents
        # and PS points to use to evaluate them. 
        hike = self.mapper.walk(PS_point, counterterm)
        # The structure of this object output should reflect each nesting level, that is:
        # hike = [ 
        #  nesting level 1. (furthest away from ME) -->   [(PS1.1, current1.1, jac1.1, kin_var1.1), (PS1.2, current1.2, jac1.1, kin_var1.2), ...],
        #  nesting level 2.                         -->   [(PS2.1, current2.1, jac2.1, kin_var2.1), (PS2.2, current2.2, jac2.2, kin_var2.2), ...],
        #  ...
        #  last nesting level (closest to ME)       -->   [(PSlast.1, currentlast.1, jaclast.1, kin_varlast.1), (PSlast.2, currentlast.2, jaclast.2, kin_varlast.2), ...],
        #  ME level                                 -->   [(BornPS, ReducedProcessInstance, None)]
        #        ]

        # Then the above "hike" can be used to evaluate the currents first and the ME last.
        # Note that the code below can become more complicated when needing to track helicities, but let's forget this for now.
        weight = 1.0
        assert((hel_config is None))
        for level, stroll in enumerate(hike[:-2]):
            for (PS_point_for_current, current, jacobian, kinematic_variables) in stroll:
                current_evaluation, all_current_results = self.all_MEAccessors(
                    current, PS_point_for_current, hel_config=None, kinematic_variables=kinematic_variables)
                # Make sure no spin- or color-correlations are demanded by the current at this stage
                assert(current_evaluation['spin_correlations']==[None,])
                assert(current_evaluation['color_correlations']==[None,])
                assert(current_evaluation['values'].keys()==[(0,0),])
                # WARNING:: this can only work for local 4D subtraction counterterms! 
                # For the integrated ones it is very likely that we cannot use a nested structure, and
                # there will be only one level anyway so that there is not need of fancy combination
                # of Laurent series.
                weight *= jacobian*current_evaluation['values'][(0,0)]['finite']
         
        # all_necesary_ME_calls is a list of tuples of the following form:
        #  (spin_correlator, color_correlator, weight)
        all_necessary_ME_calls = [(None, None, weight)]
        
        # Specify here how to combine one set of correlators with another.
        def combine_correlators(correlators_A, correlators_B):
            combined_correlator = [None, None, correlators_A[2]*correlators_B[2]]
            # Trivial combination of correlators first
            for i in range(2):
                if (correlators_A[i] is None) and (correlators_B[i] is None):
                    combined_correlator[i] = None
                elif (correlators_A[i] is None):
                    combined_correlator[i] = correlators_B[i]
                elif (correlators_B[i] is None):
                    combined_correlator[i] = correlators_A[i]
            # Non-trivial combinations now
            # Spin-correlators
            if (not correlators_A[0] is None) and (not correlators_B[0] is None):
                combined_correlator[0] = tuple(list(correlators_A[0])+list(correlators_B[0]))
            # Color-correlators
            if (not correlators_A[1] is None) and (not correlators_B[1] is None):
                # We don't know what to do yet, because we haven't agreed on how to organize
                # NNLO color correlations yet.
                raise NotImplementedError
            return combined_correlator

        # The next-to-last layer needs to be treated specifically since it must track the color- and spin-correlations        
        for (PS_point_for_current, current, jacobian, kinematic_variables) in hike[-1]:
            current_evaluation, all_current_results = self.all_MEAccessors(
                current, PS_point_for_current, hel_config=None, kinematic_variables=kinematic_variables)
            new_all_necessary_ME_calls = []
            # Now loop over all spin- and color- correlators required for this current
            # and update the necessary calls to the ME
            new_all_necessary_ME_calls = []
            for ((spin_index, color_index), current_wgt) in current_evaluation['values'].items():
                # Now combine the correlators necessary for this current, with those already
                # specified in 'all_necessary_ME_calls'
                current_correlator = ( current_evaluation['spin_correlations'][spin_index],
                                       current_evaluation['color_correlations'][color_index],
                                       jacobian*current_wgt['finite'] )
                for ME_call in all_necessary_ME_calls:
                    new_all_necessary_ME_calls.append( combine_correlators(ME_call,current_correlator) )
            # Update the list of necessary ME calls
            all_necessary_ME_calls = new_all_necessary_ME_calls

        # Finally the next layer contains the ME so it should of course be special
        assert(len(hike[-1])==1)
        ME_PS, ME_process = hike[-1]

        final_weight = 0.0        
        for (spin_correlators, color_correlators, current_weight) in all_necessary_ME_calls:
            ME_evaluation, all_ME_results = all_accessors(my_process, PS_point, alpha_s, mu_r)
            self.all_MEAccessors(
               ME_process, ME_PS, alpha_s, mu_r,
               # Let's worry about the squared order laters, we will probably directly fish
               # them out from the ME_process, since they should be set to a unique combination in
               # this case.
               squared_orders    = None,
               color_correlation = color_correlators,
               spin_correlation  = spin_correlators, 
               hel_config        = hel_config 
            )
            # Again, for the integrated subtraction counterterms, some care will be needed here
            # for the real-virtual, depending on how we want to combine the two Laurent series.
            final_weight += current_weight*ME_evaluation['finite']
        
        # Returns the corresponding weight and the mapped PS_point.
        # Also returns the mapped_process (for calling the observables), which
        # is typically simply a reference to counterterm.current which is an instance of Process.
        return final_weight, ME_PS, ME_process

    def sigma(self, PS_point, process, flavors, flavor_wgt, mu_r, mu_f1, mu_f2, *args, **opts):
        
        for ct in self.counterterms[processkey]:
            wgt += self.evaluate_counterterm(counterterm)
        
        return super(ME7Integrand_R, self).sigma(PS_point, process, flavors, flavor_wgt, mu_r, mu_f1, mu_f2, *args, **opts)

    def find_counterterms_matching_limit_type_with_regexp(self, counterterms, limit_type=None):
        """ Find all mappings that match a particular limit_type given in argument (takes a random one if left to None)."""

        # First select only the counterterms which are not pure matrix elements (i.e. they have singular structures).
        selected_counterterms = [ct for ct in counterterms if ct.is_singular()]

        returned_counterterms = []
        if not limit_type:
            returned_counterterms.append(random.choice(selected_counterterms))
        else:
            for counterterm in selected_counterterms:
                if re.match(limit_type,counterterm.get_singular_structure_string(
                                                                print_n=True, print_pdg=False, print_state=False)):
                    returned_counterterms.append(counterterm)
        
        return returned_counterterms

    def is_part_of_process_selection(self, process_list, selection=None):
        """ Checks wether any of the specified processes in the process_list provided matches the user's process
        selection. If not provided, returns True by default. 'selection' is a dictionary with the format:
           {'in_pdgs'  : ( (in_pgs1), (in_pdgs2), ...)
            'out_pdgs' : ( (out_pgs1), (out_pdgs2), ...) 
            'n_loops'  : n_loops }"""

        for process in process_list:
            is_a_match = True
            for i, in_pdg in enumerate(process.get_initial_ids()):
                if selection['in_pdgs'] is None:
                    break
                if i >= len(selection['in_pdgs']):
                    is_a_match = False
                    break
                if in_pdg not in selection['in_pdgs'][i]:
                    is_a_match = False
                    break
            if not is_a_match:
                continue
            for i, out_pdg in enumerate(process.get_final_ids_after_decay()):
                if selection['out_pdgs'] is None:
                    break
                if i >= len(selection['out_pdgs']):
                    is_a_match = False
                    break
                if out_pdg not in selection['out_pdgs'][i]:
                    is_a_match = False
                    break            
            if not is_a_match:
                continue
            if (not selection['n_loops'] is None) and process.get('n_loops') != selection['n_loops']:
                continue
            return True

        return False

    def test_IR_limits(self, test_options):
        """ Tests that 4D local subtraction terms tend to the corresponding real-emission matrix elements."""
        
        if test_options['seed']:
            random.seed(test_options['seed'])
        
        # First generate an underlying Born
        # Specifying None forces to use unformly random generating variables.
        a_real_emission_PS_point, _, _, _ = self.phase_space_generator.get_PS_point(None)

        for process_key, (defining_process, mapped_processes) in self.processes_map.items():
            # Make sure that the selected process satisfies the selected process
            if not self.is_part_of_process_selection([defining_process,]+mapped_processes, 
                                                                            selection = test_options['process']):
                continue
            
            # Here we use correction_order to select CT subset
            counterterms_to_consider = [ ct for ct in self.counterterms[process_key] if 
                        ct.count_unresolved() <= test_options['correction_order'].count('N') ]
            
            # Here we use limit_type to select the mapper to use for approaching the limit (
            # it is clear that all CT will still use their own mapper to retrieve the PS point
            # and variables to call the currents and reduced processes).
            selected_counterterms = self.find_counterterms_matching_limit_type_with_regexp(
                                                                    counterterms_to_consider, test_options['limit_type'])
            
            #misc.sprint(defining_process.nice_string())
            #misc.sprint('\n'+'\n'.join( ct.get_singular_structure_string() for ct in selected_counterterms ))
            
            # Now loop over all mappings to consider
            for limit_specifier_counterterm in selected_counterterms:
                # First identified the reduced PS point from which we can evolve to larger multiplicity
                # while becoming progressively closer to the IR limit.
                a_born_PS_point, starting_jacobian, starting_variables = self.mapper.walk_to_lower_multiplicity(
                                        a_real_emission_PS_point, limit_specifier_counterterm, kinematic_variables=True)
                
                # Now progressively approach the limit
                evaluations = {}
                # l is the scaling variable
                n_steps = test_options['n_steps']
                min_value = test_options['min_scaling_variable']
                acceptance_threshold = test_options['acceptance_threshold']

                for scaling_parameter in range(1,n_steps+1):
                    # Use equally spaced steps on a log scale
                    scaling_parameter = 10.0**(-((float(scaling_parameter)/n_steps)*abs(math.log10(min_value))))
                    scaled_real_PS_point, jacobian = self.mapper.approach_limit(
                                a_born_PS_point, limit_specifier_counterterm, starting_variables, scaling_parameter)
                    
                    # Evaluate  real ME
                    assert(len([ct for ct in self.counterterms[process_key] if not ct.is_singular()])==0)
                    non_singular_ME = [ct for ct in self.counterterms[process_key] if not ct.is_singular()][0]

                    # Need smarter way to evaluate the ME
                    ME_evaluation, _, _ = self.evaluate_counterterm(non_singular_ME, scaled_real_PS_point, hel_config=None)
                    
                    # Approximated real ME (aka. local 4d subtraction counterterm)
                    summed_counterterm_weight = 0.0
                    for counterterm in counterterms_to_consider:
                        if counterterm.is_singular():
                            continue
                        ct_weight, _, _ = self.evaluate_counterterm(counterterm, scaled_real_PS_point, hel_config=None)
                        summed_counterterm_weight += ct_weight
                    
                    # Add evaluations to the list so as to study how the approximated reals converge towards the real
                    evaluations[scaling_parameter]= {
                         'non_singular_ME'      : ME_evaluation, 
                         'approximated_ME'      : summed_counterterm_weight,
                         'limit_specifier'      : limit_specifier_counterterm,
                         'defining_process'     : defining_process
                        }
                
                all_evaluations[(process_key, counterterm.get_singular_structure_string(
                                                        print_n=True, print_pdg=False, print_state=False))] = evaluations

        # Now produce a nice matplotlib of the evaluations and assess whether this test passed or not.
        return self.analyze_IR_limits_test(all_evaluations, acceptance_threshold)

    def analyze_IR_limits_test(self, all_evaluations, acceptance_threshold):
        """ Analyze the results of the test_IR_limits command. """
        #TODO
        pass


class ME7Integrand_RR(ME7Integrand_R):
    """ ME7Integrand for the computation of a double real-emission type of contribution."""
    def sigma(self, PS_point, process, flavors, flavor_wgt, mu_r, mu_f1, mu_f2, *args, **opts):
        return super(ME7Integrand_RR, self).sigma(PS_point, process, flavors, flavor_wgt, mu_r, mu_f1, mu_f2, *args, **opts)

# Integrand classes map is defined here as module variables. This map can be overwritten
# by the interface when using a PLUGIN system where the user can define his own Integrand.
# Notice that this must be placed after all the Integrand daughter classes in this module have been declared.
ME7Integrand_classes_map = {'Born': ME7Integrand_B,
                            'LoopInduced_Born': ME7Integrand_LIB,
                            'Virtual': ME7Integrand_V,
                            'SingleReals': ME7Integrand_R,
                            'DoubleReals': ME7Integrand_RR,
                            'Unknown': None}