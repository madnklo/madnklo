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
import itertools

import stat
import subprocess
import sys
import time
from datetime import datetime
import tarfile
import StringIO
import shutil
import copy
import json
from collections import OrderedDict
from distutils.util import strtobool

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


import madgraph.core.base_objects as base_objects
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
import madgraph.integrator.ME7_integrands as ME7_integrands
#    import madgraph.various.histograms as histograms  # imported later to not slow down the loading of the code
import models.check_param_card as check_param_card
import models.model_reader as model_reader
import models.import_ufo as import_ufo

from madgraph.iolibs.files import ln
from madgraph import InvalidCmd, MadGraph5Error, MG5DIR, ReadWrite, MPI_RANK, MPI_SIZE, MPI_ACTIVE



class MadEvent7Error(Exception):
    pass

class ZeroResult(MadEvent7Error):
    pass

#===============================================================================
# Utility
#===============================================================================

# To be used as a with statement
class ME7RunEnvironment(misc.Silence, misc.MuteLogger):
    """ Sets up an environment optimized for the various runs of ME7 interface:
      > The ProcessKey optimization is temporarily activated
      > The stdout/stderr is temporarily redirected to \dev\null"""

    def __init__(self, 
                 silence=True, 
                 accessor_optimization=True, 
                 loggers = [('madgraph',logging.INFO),
                            ('contributions',logging.INFO),
                            ('madevent7',logging.INFO),
                            ('madevent7.stderr',logging.INFO)]
                ):
        self.silence = silence
        if self.silence:
            misc.Silence.__init__(self)
        self.accessor_optimization = accessor_optimization
        # Allow for initialization with just a logger level value
        if isinstance(loggers, int):
            loggers = [('madgraph',loggers),
                       ('contributions',loggers),
                       ('madevent7',loggers),
                       ('madevent7.stderr',loggers)]
        self.logger_names = [pair[0] for pair in loggers]
        self.logger_levels = [pair[1] for pair in loggers]
        if len(self.logger_names)>0:
            misc.MuteLogger.__init__(self, self.logger_names, self.logger_levels)
        
    def __enter__(self):
        if len(self.logger_names)>0:
            misc.MuteLogger.__enter__(self)
        if self.silence:
            misc.Silence.__enter__(self)
        from madgraph.core.accessors import activate_cache
        if self.accessor_optimization:
            activate_cache()

    def __exit__(self, *args):
        if len(self.logger_names)>0:
            misc.MuteLogger.__exit__(self, *args)
        if self.silence:
            misc.Silence.__exit__(self, *args)
        from madgraph.core.accessors import deactivate_cache
        if self.accessor_optimization:
            deactivate_cache()

# To be used as a decorator
def wrap_with_ME7RunEnvironment(*decorator_args, **decorator_opts):
    """ Decorate a function and automatically enclose it within a with statement
    providing the ME7Run environment."""
    def add_environment_to_function(f):
        
        def modified_function(*args, **opt):
            with ME7RunEnvironment(*decorator_args,**decorator_opts):
                return f(*args, **opt)
        return modified_function
    
    return add_environment_to_function

def mute_logger(names=['madgraph','ALOHA','cmdprint','madevent'], levels=[50,50,50,50]):
    """change the logger level and restore those at their initial value at the
    end of the function decorated."""
    def control_logger(f):
        def restore_old_levels(names, levels):
            for name, level in zip(names, levels):
                log_module = logging.getLogger(name)
                log_module.setLevel(level)            
        
        def f_with_no_logger(self, *args, **opt):
            old_levels = []
            for name, level in zip(names, levels):
                log_module = logging.getLogger(name)
                old_levels.append(log_module.level)
                log_module.setLevel(level)
            try:
                out = f(self, *args, **opt)
                restore_old_levels(names, old_levels)
                return out
            except:
                restore_old_levels(names, old_levels)
                raise
            
        return f_with_no_logger
    return control_logger

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
    
    def parse_test_IR_options(self, args, mode='limits'):
        """ Parses the options specified to the test_limit command."""
        
        assert(mode in ['limits', 'poles'])

        # None means unspecified, therefore considering all types.
        testlimits_options = {
            'correction_order'        : None,
            'limits'                  : [None,],
            'counterterms'            : None,
            'walker'                  : 'LorentzNLO',
            'process'                 : {'in_pdgs'  : None,
                                         'out_pdgs' : None,
                                         'n_loops'  : None},
            'seed'                    : None,
            'n_steps'                 : 30,
            'max_scaling_variable'    : 1.,
            'min_scaling_variable'    : 1.0e-6,
            'acceptance_threshold'    : 1.0e-4,
            'apply_higher_multiplicity_cuts' : True,
            'apply_lower_multiplicity_cuts'  : True,
            'show_plots'              : True,
            'save_plots'              : False,
            'save_results_to_path'    : None,
            'plots_suffix'            : None,
            'ignore_flavors'          : False,
            'set_PDFs_to_unity'       : True,
            'boost_back_to_com'       : True,
            'epsilon_expansion_term'  : 'sum_all',
            'selected_sectors'        : None,
            # Here we store a list of lambda function to apply as filters
            # to the integrand we must consider
            'integrands'              : [lambda integrand: True]
        }

        if mode=='poles':
            # Remove some options for the pole
            del testlimits_options['limits']
            del testlimits_options['walker']
            del testlimits_options['n_steps']
            del testlimits_options['min_scaling_variable']
            del testlimits_options['max_scaling_variable']
            del testlimits_options['show_plots']
            del testlimits_options['save_plots']
            del testlimits_options['plots_suffix']
            del testlimits_options['boost_back_to_com']
            # We must be much tighter in the check of the relative difference of the pole
            # residues.
            testlimits_options['acceptance_threshold'] = 1.0e-13
        elif mode=='limits':
            pass   
 
        # Group arguments in between the '--' specifiers.
        # For example, it will group '--process=p' 'p' '>' 'd' 'd~' 'z'
        # into '--process=p p > d d~ z'.
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
            elif key in ['--apply_cuts', '--ac']:
                if value is None or value.lower() in ['t','true']:
                    testlimits_options['apply_higher_multiplicity_cuts'] = True
                    testlimits_options['apply_lower_multiplicity_cuts'] = True
                elif value.lower() in ['f','false']:
                    testlimits_options['apply_higher_multiplicity_cuts'] = False
                    testlimits_options['apply_lower_multiplicity_cuts'] = False 
                elif value.lower()=='no_lower':
                    testlimits_options['apply_lower_multiplicity_cuts'] = False
                elif value.lower()=='no_higher':
                    testlimits_options['apply_higher_multiplicity_cuts'] = False                             
                elif value.lower()=='lower':
                    testlimits_options['apply_lower_multiplicity_cuts'] = True
                elif value.lower()=='higher':
                    testlimits_options['apply_higher_multiplicity_cuts'] = True
                else:
                    raise InvalidCmd("Value '%s' for option '%s' not recognized."%(value,key))
            elif key in ['--n_steps'] and mode=='limits':
                try:
                    testlimits_options[key[2:]] = int(value)
                    if int(value)<2:
                        raise ValueError  
                except ValueError:
                    raise InvalidCmd("'%s' is not a valid integer for option '%s'"%(value, key))
            elif key == '--save_results_to_path':
                if value == 'None':
                    testlimits_options['save_results_to_path'] = None
                else:
                    if os.path.isabs(value):
                        path = value
                    else:
                        path = pjoin(self.me_dir, value)
                    if os.path.isfile(path):
                        logger.warning("File path '%s' for "%path+
                             "saving test results already exists and will be overwritten.")
                        os.remove(path)
                    testlimits_options['save_results_to_path'] = path
            elif key in ['--counterterms']:
                if not isinstance(value, str):
                    raise InvalidCmd("'%s' is not a valid option for '%s'"%(value, key))
                testlimits_options[key[2:]] = value
            elif key in ['--epsilon_expansion_term','--eet']:
                if value.lower() not in ['sum_all','finite'] and re.match('^eps\^(-?\d*)$', value.lower().strip()) is None:
                    raise InvalidCmd("'%s' is not a valid option for '%s'. It must be 'finite', 'sum_all' or of the form 'eps^<i>'."%(value, key))
                testlimits_options['epsilon_expansion_term'] = value
            elif key in ['--acceptance_threshold']:
                try:
                    testlimits_options[key[2:]] = float(value)                  
                except ValueError:
                    raise InvalidCmd("'%s' is not a valid float for option '%s'"%(value, key))
            elif key in ['--min_scaling_variable','--ms'] and mode=='limits':
                try:
                    testlimits_options['min_scaling_variable'] = float(value)                  
                except ValueError:
                    raise InvalidCmd("'%s' is not a valid float for option '%s'"%(value, key))
            elif key in ['--max_scaling_variable'] and mode == 'limits':
                try:
                    testlimits_options['max_scaling_variable'] = float(value)
                except ValueError:
                    raise InvalidCmd("'%s' is not a valid float for option '%s'"%(value, key))
            elif key in ['--walker'] and mode=='limits':
                try:
                    testlimits_options['walker'] = value
                except KeyError:
                    raise InvalidCmd("'%s' is not a valid %s" % (value, key[2:]))
            elif key in ['--limits','--l'] and mode=='limits':
                if not isinstance(value, str):
                    raise InvalidCmd("'%s' is not a valid option for '%s'"%(value, key))
                # Check if a list of defining limits is specified
                try:
                    evaluated_list = eval(value)
                    if not isinstance(evaluated_list, list):
                        raise BaseException
                    testlimits_options['limits'] = evaluated_list
                except:
                    testlimits_options['limits'] = [value,]

            elif key in ['--selected_sector', '--sector', '--ss']:
                eval_value = eval(value)
                if testlimits_options['selected_sectors'] is None:
                    if isinstance(eval_value, list):
                        testlimits_options['selected_sectors'] = eval_value
                    else:
                        testlimits_options['selected_sectors'] = [ eval_value,]
                else:
                    if isinstance(eval_value, list):
                        testlimits_options['selected_sectors'].extend(eval_value)
                    else:
                        testlimits_options['selected_sectors'].append(eval_value)

            elif key in [
                '--show_plots','--set_PDFs_to_unity', '--save_plots',
                '--boost_back_to_com', '--ignore_flavors' ]:
                try:
                    if value is None:
                        testlimits_options[key[2:]] = True
                    else:
                        testlimits_options[key[2:]] = strtobool(value)
                except ValueError:
                    raise InvalidCmd("Cannot set '%s' option to '%s'."%(key, value))
            elif key == '--seed':
                try:
                    testlimits_options['seed'] = int(value)
                except ValueError:
                    raise InvalidCmd("Cannot set '%s' option to '%s'."%(key, value))
            elif key == '--plots_suffix':
                testlimits_options['plots_suffix'] = value
            elif key in ['--integrands', '--itg']:
                testlimits_options['integrands'].extend(self.get_integrand_filters(value, 'select'))
            elif key in ['--veto_integrands', '--veto_itg']:
                testlimits_options['integrands'].extend(self.get_integrand_filters(value, 'reject'))
            elif key=='--process':
                if value.lower()=='all':
                    testlimits_options['process']['in_pdgs'] = None
                    testlimits_options['process']['out_pdgs'] = None
                else:
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
                                raise InvalidCmd("Particle '%s' not recognized in current model."%part)
                        final_pdgs.append(tuple(multi_part))
                        
                    testlimits_options['process']['in_pdgs'] = tuple(initial_pdgs)
                    testlimits_options['process']['out_pdgs'] = tuple(final_pdgs)
            elif key == '--n_loops':
                if value.lower()=='all':
                    testlimits_options['process'][key[2:]] = None
                else:
                    try:
                        testlimits_options['process'][key[2:]] = int(value)
                        if int(value)<0:
                            raise ValueError
                    except ValueError:
                        raise InvalidCmd("'%s' is not a valid integer for option '%s'"%(value, key))

            else:
                raise InvalidCmd("Option '%s' for the test_IR_%s command not recognized."%(key,mode))        
        
        return new_args, testlimits_options

    def parse_set_integrator_options(self, args):
        """ Parses the argument of the set_integrator command."""

        # First combine all value of the options (starting with '--') separated by a space
        opt_args = []
        new_args = []
        for arg in args:
            if arg.startswith('--'):
                opt_args.append(arg)
            elif len(opt_args) > 0:
                opt_args[-1] += ' %s' % arg
            else:
                new_args.append(arg)

        set_integrator_options = {
                                  # Vegas3 options
                                  'n_points'            : 1000,
                                  'n_points_survey'     : 1000,
                                  'n_points_refine'     : 5000,
                                  'n_iterations'        : 10,
                                  'n_iterations_survey' : 10,
                                  'n_iterations_refine' : 10,                          
                                  'verbosity'           : 1,
                                  'batch_size'          : 1000,
                                  'seed'                : None,
                                  'save_grids'          : None,
                                  'load_grids'          : None,
                                  # Cuba@Vegas options
                                  'max_eval'            : int(1e10),
                                  'max_eval_survey'     : 10000,
                                  'n_start'             : 1000,
                                  'n_start_survey'      : 1000,
                                  'n_increase'          : 1000,
                                  'n_increase_survey'   : 1000,
                                  'target_accuracy'     : 1.0e-3,
                                  'target_accuracy_survey' : 5.0e-2
                                  }

        for arg in opt_args:
            try:
                key, value = arg.split('=',1)
            except ValueError:
                key = arg
                value = None
            
            if key == '--integrator':
                if value.upper() not in self._integrators:
                    raise InvalidCmd("Selected integrator '%s' not recognized."%value)
                set_integrator_options['integrator'] = value.upper()
            elif key in ['--n_points']:
                set_integrator_options[key[2:]] = int(value)
                set_integrator_options[key[2:]+'_survey'] = int(value)
                set_integrator_options[key[2:]+'_refine'] = 5*int(value)
            elif key in ['--n_iterations']:
                set_integrator_options[key[2:]] = int(value)
                set_integrator_options[key[2:]+'_survey'] = int(value)
                set_integrator_options[key[2:]+'_refine'] = int(value)
            elif key in ['--n_points_survey', '--n_iterations_survey',
                         '--n_points_refine', '--n_iterations_refine']:
                set_integrator_options[key[2:]] = int(value)
            elif key=='--verbosity':
                modes = {'none':0, 'integrator':1, 'all':2}
                set_integrator_options[key[2:]] = modes[value.lower()]
            elif key in ['--seed','--max_eval','--max_eval_survey','--n_start',
                         '--n_start_survey','--n_increase','--n_increase_survey']:
                try:
                    set_integrator_options[key[2:]] = int(value)
                except ValueError:
                    raise InvalidCmd("Cannot set '%s' option to '%s'."%(key, value))
            elif key in ['--target_accuracy','--target_accuracy_survey']:
                try:
                    set_integrator_options[key[2:]] = float(value)
                except ValueError:
                    raise InvalidCmd("Cannot set '%s' option to '%s'."%(key, value))
            elif key in ['--save_grids', '--load_grids']:
                set_integrator_options[key[2:]] = value
            elif key in ['--batch_size','--bs']:
                try:
                    set_integrator_options['batch_size'] = int(value)
                except ValueError:
                    raise InvalidCmd("Value '%s' not valid integer for option '%s'."%(value,key))                
            else:
                # We have not coded here all the options the integrator supports, so just
                # blindly assign them for now
                set_integrator_options[key[2:]] = eval(value)

        return new_args, set_integrator_options

    def parse_launch(self, args):
        """ Parses the argument of the launch command."""

        # First combine all value of the options (starting with '--') separated by a space
        opt_args = []
        new_args = []
        for arg in args:
            if arg.startswith('--'):
                opt_args.append(arg)
            elif len(opt_args) > 0:
                opt_args[-1] += ' %s' % arg
            else:
                new_args.append(arg)

        launch_options = {'integrator'          : 'VEGAS3',
                          'verbosity'           : 1,
                          'refresh_filters'     : 'auto',
                          'compile'             : 'auto',
                          'seed'                : None,
                          # Here we store a list of lambda function to apply as filters
                          # to the integrand we must consider
                          'integrands'          : [lambda integrand: True],
                          'run_name'            : ''
                          }
        
        for arg in opt_args:
            try:
                key, value = arg.split('=',1)
            except ValueError:
                key = arg
                value = None
            
            if key == '--integrator':
                if value.upper() not in self._integrators:
                    raise InvalidCmd("Selected integrator '%s' not recognized."%value)
                launch_options['integrator'] = value.upper()
            elif key=='--verbosity':
                modes = {'none':0, 'integrator':1, 'all':2}
                launch_options[key[2:]] = modes[value.lower()]
            elif key == '--seed':
                try:
                    launch_options['seed'] = int(value)
                except ValueError:
                    raise InvalidCmd("Cannot set '%s' option to '%s'."%(key, value))
            elif key in ['--refresh_filters','--compile']:
                available_modes = ['auto','never','always']
                if value is None:
                    mode = 'always'
                else:
                    mode = str(value).lower()
                if mode not in available_modes:
                    raise InvalidCmd("Value '%s' not valid for option '%s'."%(value,key))
                launch_options[key[2:]] = mode

            elif key in ['--integrands', '--itg']:
                launch_options['integrands'].extend(self.get_integrand_filters(value, 'select'))
            
            elif key in ['--veto_integrands', '--veto_itg']:
                launch_options['integrands'].extend(self.get_integrand_filters(value, 'reject'))

            elif key=='--run_name':
                launch_options['run_name'] = value
            else:
                raise InvalidCmd("Option '%s' for the launch command not recognized."%key)

        return new_args, launch_options

    def parse_display_integrands(self, args):
        """ Parses the argument of the "display contributions" command."""

        display_options = {'format':0,
                          # Here we store a list of lambda function to apply as filters
                          # to the ingegrand we must consider
                          'integrands': [lambda integrand: True]}        
        
        for arg in args:
            try:
                key, value = arg.split('=',1)
            except ValueError:
                key = arg
                value = None
            
            if key in ['--format']:
                try:
                    display_options['format'] = int(value)
                    if display_options['format']<0:
                        raise ValueError
                except ValueError:
                    raise InvalidCmd("Value '%s' invalid for option 'format'."%value)
            elif key in ['--integrands', '--itg']:
                display_options['integrands'].extend(self.get_integrand_filters(value, 'select'))
            
            elif key in ['--veto_integrands', '--veto_itg']:
                display_options['integrands'].extend(self.get_integrand_filters(value, 'reject'))
            
            else:
                raise InvalidCmd("Option '%s' for the display integrands command not recognized."%key)

        return display_options

    def get_integrand_filters(self, filters, mode):
        """ Given user defined filters (a string), returns a lambda function on
        a ME7 integrand that applies it. mode can either be 'select' or 'reject',
        depending on whether one wants to keep or reject those integrands."""
        
        assert(mode in ['select','reject'])
        
        if filters is None:
            InvalidCmd("Option '--integrands' must specify values with '='.")
        
        expansion_orders = []
        integrand_types  = []
        short_names      = []
        
        integrand_type_short_cuts = dict(
            (k, eval('ME7_integrands.ME7Integrand_%s'%k, ))
            for k in ['R','RV','RR','RRR','V','B'] )
        
        for filter in filters.split(','):
            f = filter.strip()
            if f in integrand_type_short_cuts:
                integrand_types.append(integrand_type_short_cuts[f])
            elif f in ['LO','NLO','NNLO','NNNLO']:
                expansion_orders.append(f)
            else:
                try:
                    integrand_type = eval(f)
                    if isinstance(integrand_type, type):
                        integrand_types.append(integrand_type)
                    else:
                        raise BaseException
                except:
                    short_names.append(f)
        
        filter_functions = []
        if mode == 'select':
            # First the ordering filter if specified
            if expansion_orders:
                filter_functions.append(lambda integrand: 
                    integrand.contribution_definition.correction_order in expansion_orders)
            # The the expansion orders
            if integrand_types:
                filter_functions.append(lambda integrand: (type(integrand) in integrand_types) )
            if short_names:
                filter_functions.append(lambda integrand: (
                          integrand.contribution_definition.short_name() in short_names) )                
        else:
            # First the ordering filter if specified
            if expansion_orders:
                filter_functions.append(lambda integrand: 
                    not (integrand.contribution_definition.correction_order in expansion_orders) )
            # Then the expansion orders
            if integrand_types:
                filter_functions.append(lambda integrand: (type(integrand) not in integrand_types) )
            if short_names:
                filter_functions.append(lambda integrand: (
                        integrand.contribution_definition.short_name() not in short_names) )
        
        return filter_functions

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
    integrator_verbosity = 1 if logger.level > logging.DEBUG else 2
    _integrators = {
       'NAIVE' : (integrators.SimpleMonteCarloIntegrator, 
                  { 'n_iterations'            : 10,
                    'n_points_per_iterations' : 100,
                    'accuracy_target'         : None,
                    'verbosity'               : integrator_verbosity  } ),
    
       'VEGAS3' : (vegas3_integrator.Vegas3Integrator,
                   { 'n_iterations_survey'     : 10,
                     'n_points_survey'         : 1000,
                     'n_iterations_refine'     : 10,
                     'n_points_refine'         : 2000,
                     'verbosity'               : integrator_verbosity  } ),
    
       'VEGAS' : (pyCubaIntegrator.pyCubaIntegrator, 
                  { 'algorithm' : 'Vegas', 
                    'verbosity' : integrator_verbosity,
                    'seed'      : 0,
                    'target_accuracy' : 1.0e-3,
                    'target_accuracy_survey' : 1.0e-3,
                    'n_start'           : 1000,
                    'n_start_survey'    : 1000,
                    'n_increase'        : 500,
                    'n_increase_survey' : 500,
                    'n_batch'           : 1000,
                    'max_eval'          : int(1e10),
                    'max_eval_survey'   : 10000,
                    'min_eval'          : 0}),
       
       'SUAVE'   : (pyCubaIntegrator.pyCubaIntegrator, 
                    { 'algorithm' :'Suave', 
                      'verbosity' : integrator_verbosity,
                      'target_accuracy' : 1.0e-3,
                      'max_eval'  : int(1e10),
                      'min_eval'  : 0 } ),
      
       'DIVONNE' : (pyCubaIntegrator.pyCubaIntegrator, 
                    { 'algorithm' : 'Divonne', 
                      'verbosity': integrator_verbosity,
                      'target_accuracy' : 1.0e-5,
                      'max_eval'  : int(1e10),
                      'min_eval'  : 0 } ),
    
       'CUHRE'   : (pyCubaIntegrator.pyCubaIntegrator, 
                    { 'algorithm' : 'Cuhre',
                      'verbosity' : integrator_verbosity,
                      'target_accuracy' : 1.0e-3,
                      'max_eval'  : int(1e10),
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
        self.all_integrands = ME7_integrands.ME7IntegrandList()
        self.model = None
        self.run_card = None
        # Overwrite the above properties from the bootstrap using the database written at the output stage
        self.bootstrap_ME7()
        
        # Instance of the current integrator
        self.integrator = None

    def check_output_type(self, path):
        """ Check that the output path is a valid Madevent 7directory """        
        return os.path.isfile(os.path.join(path,'MadEvent7.db'))

    def check_already_running(self):
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
            ME7_dump['all_MEAccessors'], root_path=self.me_dir, model=self.model)

        # Pass on the process directory
        self.options['me_dir'] = self.me_dir

        self.all_integrands = ME7_integrands.ME7IntegrandList([
            integrand_dump['class'].initialize_from_dump(
                integrand_dump,
                self.model, self.run_card, self.all_MEAccessors, self.options)
            for integrand_dump in ME7_dump['all_integrands'] ])
        
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
        # Check if run_card values are supported.
        self.run_card.check_validity()

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

    def do_display(self, line, *args, **opts):
        """ General display command, the first argument of the line should be what must
        be displayed."""

        try:
            object_to_display = self.split_arg(line)[0]
            display_function = getattr(self,'display_%s'%object_to_display)
        except IndexError, AttributeError:
            return super(MadEvent7Cmd, self).do_display(line, *args, **opts)
        
        display_function(line, *args, **opts)

    def display_integrands(self, line, *args, **opts):
        """ Displays integrands."""
        
        args = self.split_arg(line)
        display_options = self.parse_display_integrands(args[1:])
        
        integrands_to_consider = ME7_integrands.ME7IntegrandList([ itg for itg in self.all_integrands if
                           all(filter(itg) for filter in display_options['integrands']) ])

        logger.info('\n'+integrands_to_consider.nice_string(format=display_options['format']))

    def do_set_integrator_options(self, line):
        """ Command allowing to specify options for a specific integrator."""
        
        args = self.split_arg(line)
        new_args, integrator_options = self.parse_set_integrator_options(args)

        if len(new_args)==0:
            raise InvalidCmd('An integrator in %s must be specified as the first '%str(self._integrators)+
                                         'argument in the command set_integrator_option.')
        elif new_args[0] not in self._integrators:
            raise InvalidCmd("The specified integrator '%s' is not in the list of supported ones (%s)."%(
                                                      new_args[0], str(self._integrators.keys())))            

        self._integrators[new_args[0]][1].update(integrator_options)

    def do_launch(self, line, *args, **opt):
        """Main command, starts the cross-section computation. Very basic setup for now.
        We will eventually want to have all of these meta-data controllable via user commands
        in the interface.
        This is super naive and only for illustrative purposes for now."""
        
        args = self.split_arg(line)
        new_args, launch_options = self.parse_launch(args)

        # In principle we want to start by recompiling the process output
        # to make sure that everything is up to date.
        self.synchronize(**launch_options)

        # Setup parallelization
        self.configure_run_mode(self.options['run_mode'])

        # Re-initialize the integrator
        integrator_name = launch_options['integrator']
        integrator_options = self._integrators[integrator_name][1]

        integrator_options['verbosity'] = launch_options['verbosity']
        integrator_options['cluster'] = self.cluster

        if integrator_name=='VEGAS3':
            integrator_options['parallelization'] = self.cluster

        integrands_to_consider = ME7_integrands.ME7IntegrandList([
            itg for itg in self.all_integrands
            if all(filter(itg) for filter in launch_options['integrands']) ])
        self.integrator = self._integrators[integrator_name][0](
            integrands_to_consider, **integrator_options)
        
        if len(set([len(itgd.get_dimensions()) for itgd in integrands_to_consider]))>1 and integrator_name not in ['NAIVE']:
            # Skip integrators that do not support integrands with different dimensions.
            # Of course, in the end we wil not naively automatically put all integrands alltogether in the same integrator
            raise InvalidCmd("For now, whenever you have several integrands of different dimensions, only "+
                             "the NAIVE integrator is available.")

        #Setup the output folder
        if MPI_RANK==0:
            if launch_options['run_name']:
                run_output_path = pjoin(self.me_dir, 'Results', 'run_%s' % launch_options['run_name'])
            else:
                run_output_path = pjoin(self.me_dir, 'Results', 'run_%s' % self.run_card['run_tag'])
            suffix_number = 1
            suffix = '_%d' % suffix_number
            existing_results = {}

            # If a folder for this run_name exists, go through appending/creating new folder logic
            if (os.path.exists(run_output_path + suffix)):
                # Go through existing folders to find the last one
                while os.path.exists(run_output_path + suffix):
                    suffix_number += 1
                    suffix = '_%d' % suffix_number
                suffix_number -= 1
                suffix = '_%d' % suffix_number
                # Get the names of all the integrands to be integrated
                contribution_names = [integrand.contribution_definition.short_name() for integrand in self.integrator.integrands]
                # Import all results already in the data file
                try:
                    with open(pjoin(run_output_path + suffix, "cross_sections.json"), "r") as result_file:
                        existing_results = json.load(result_file)
                except (IOError,ValueError) as e:
                    logger.warn("Could not open the result file in the existing result folder %s."%(run_output_path + suffix))
                # If the importation failed (inexistant file or json problem) or if there is overlap with existing data, create a new folder
                if (not existing_results) or any([name in existing_results for name in contribution_names]):
                    logger.info("Current contribution already in the existing result folder %s."%(run_output_path + suffix))
                    suffix_number += 1
                    suffix = '_%d' % suffix_number
                    existing_results = {}
                # If the JSON datafile was correctly imported AND the integrands are not already in this data file, append there
                else:
                    logger.info("Current contribution added in the existing result folder %s."%(run_output_path + suffix))
            else:
                logger.info("This is the first run with this name")
            # Suffix is either the last existing if appending is possible or last+1 if new folder needed
            run_output_path = run_output_path + suffix
            try:
                # This fails if the folder exists already. Creation only goes through if we need it.
                os.makedirs(run_output_path)
                logger.info("Creating a new result folder at %s"%run_output_path) # Message only displayed if the new folder is created
            except OSError:
                pass

        # The integration can be quite verbose, so temporarily setting their level to 50 by default is best here
        if launch_options['verbosity'] > 1:
            logger_level = logging.DEBUG
        else:
            logger_level = logging.INFO

        logger.info("="*100)        
        logger.info('{:^100}'.format("Starting integration, lay down and enjoy..."),'$MG:color:GREEN')
        logger.info("="*100)

        # Wrap the call in a propice environment for the run
        with ME7RunEnvironment( silence = False, accessor_optimization = True, loggers = logger_level ):
            xsec, error = self.integrator.integrate()

        logger.info("="*100)
        logger.info('{:^100}'.format("Cross-section of %s@%s with integrator '%s':"%(
            os.path.basename(pjoin(self.me_dir)), '+'.join(itg.get_short_name() for itg in integrands_to_consider),
                                            self.integrator.get_name())),'$MG:color:GREEN')
        logger.info('{:^100}'.format("%.5e +/- %.2e [pb]"%(xsec, error)),'$MG:color:BLUE')
        logger.info("="*100)

        # Deal with normalizing the plotting and outputting the plots
        plot_collector = {} # We set it so that plot_collector[observable_name] is a HwU with the sum for all integrands
        # Ultimately plot_collector should be an object inside a bigger DataCollector object which allows to navigate
        # between the different results etc, which would be good to be able to sum stuff (like R+V after integration)
        for integrand in self.integrator.integrands:
            n_integrand_calls = integrand.n_observable_calls
            if n_integrand_calls <= 0:
                continue
            if not integrand.apply_observables:
                continue
            for observable in integrand.observable_list:
                
                # Assign a MC uncertainty to the plots
                observable.HwU.set_statistical_uncertainty()
                observable.HwU *= self.integrator.observables_normalization(n_integrand_calls)
                if observable.name in plot_collector:
                    plot_collector[observable.name] += observable.HwU
                else:
                    plot_collector[observable.name] = observable.HwU

        if len(plot_collector) > 0:
            # Eventually the plot collector module will handle the histograms dependencies itself.
            from madgraph.various.histograms import HwUList
            final_hwu_list = HwUList(plot_collector.values())
            HwU_plots_path = pjoin(run_output_path,'%s_plots'%self.run_card['fo_analysis'])
            logger.info("Plots from the fixed-order analysis '%s' written out to '%s.gnuplot'."%
                                             (self.run_card['fo_analysis'],HwU_plots_path))
            final_hwu_list.output(
                HwU_plots_path, 
                format='gnuplot'
            )

        logger.info('')

        if MPI_RANK==0:
            # Write the result in 'cross_sections.json' of the result directory
            self.output_run_results(
                run_output_path, xsec, error, self.integrator.integrands, existing_results)


    def do_test_IR_limits(self, line, *args, **opt):
        """This function test that local subtraction counterterms match
        with the actual matrix element in the IR limit."""
    
        args = self.split_arg(line)
        args, testlimits_options = self.parse_test_IR_options(args, mode='limits')
        
        # In principle we want to start by recompiling the process output
        # so as to make sure that everything is up to date.
        self.synchronize(**testlimits_options)
        
        if testlimits_options['correction_order'] is None:
            # If not defined, automatically assign correction_order
            # to the highest correction considered.
            testlimits_options['correction_order'] = self.mode

        n_integrands_run = 0
        for integrand in self.all_integrands:
            if not hasattr(integrand, 'test_IR_limits') or \
                not all(filter(integrand) for filter in testlimits_options['integrands']):
                    continue
            n_integrands_run += 1
            logger.debug(
                'Now testing IR limits of the following integrand:\n' +
                integrand.nice_string() )
            # Adjust the selection of sectors if necessary
            old_selector = integrand.is_sector_selected
            if testlimits_options['selected_sectors'] is not None:
                integrand.is_sector_selected = lambda defining_process, sector: \
                                            sector.leg_numbers in testlimits_options['selected_sectors']
            integrand.test_IR_limits(test_options=testlimits_options)
            integrand.is_sector_selected = old_selector

        if n_integrands_run == 0:
            logger.warning("No available integrand for function 'test_IR_limits'")

    def do_test_IR_poles(self, line, *args, **opt):
        """This function tests that the integrated counterterms match with the IR poles
        from the virtual contributions and PDF counterterms."""
    
        args = self.split_arg(line)
        args, testpoles_options = self.parse_test_IR_options(args, mode='poles')
        
        # In principle we want to start by recompiling the process output
        # so as to make sure that everything is up to date.
        self.synchronize(**testpoles_options)
        
        if testpoles_options['correction_order'] is None:
            # If not defined, automatically assign correction_order to the highest correction considered.
            testpoles_options['correction_order'] = self.mode

        n_integrands_run = 0
        for integrand in self.all_integrands:
            if not hasattr(integrand, 'test_IR_poles') or \
               not all(filter(integrand) for filter in testpoles_options['integrands']):
                continue
            n_integrands_run += 1
            logger.debug('Now testing IR poles of the following integrand:\n%s'%(integrand.nice_string()))
            # Adjust the selection of sectors if necessary
            old_selector = integrand.is_sector_selected
            if testpoles_options['selected_sectors'] is not None:
                integrand.is_sector_selected = lambda defining_process, sector: \
                                            sector.external_legs in testpoles_options['selected_sectors']
            integrand.test_IR_poles(test_options=testpoles_options)
            integrand.is_sector_selected = old_selector

        if n_integrands_run==0:
            logger.warning("No available integrand for function 'test_IR_poles'")

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

    @staticmethod
    def output_run_results(run_output_path, xsec, error, integrands, existing_results):
        """Output the run results as a json file

        Take the preexisting data in a JSON output file, add the information for the new run and update the output file with all the data. When multiple integrands are considered in one integration, we output their sum, tagged as "Contrib1+Contrib2+..." (e.g. V+R+BS).

        Arguments:
            run_output_path {[str]} -- where to save the file
            xsec {[float]} -- cross section for the set of integrands considered
            error {[float]} -- Monte Carlo uncertainty on xsec
            integrands {[ME7IntegrandList]} -- List of ME7Integrand objects that were integrated
            existing_results dict -- The JSON data for results already present in run_output_path/cross_sections.json
            """
        with open(pjoin(run_output_path, 'cross_sections.json'), 'w') as xsec_summary:
            logger.info("Saving the integration result to:")
            logger.info("%s"%(pjoin(run_output_path, 'cross_sections.json')))
            contribution_name = "+".join([integrand.contribution_definition.short_name() for integrand in integrands])
            result = {contribution_name: OrderedDict()}
            result[contribution_name]["name"] = contribution_name
            result[contribution_name]["Cross section"] = xsec
            result[contribution_name]["MC uncertainty"] = error
            result[contribution_name]["unit"] = "pb"
            orders = [integrand.contribution_definition.correction_order for integrand in integrands]
            # The overall order of the contribution is taken as the maximum of the orders considered:
            # B+V+R is NLO, RR+B is NNLO etc
            result[contribution_name]["order"] = max(orders)
            result[contribution_name]["timestamp"] = datetime.now().strftime('%Y:%m:%d:%H:%M:%S')
            result.update(existing_results)
            json.dump(result, xsec_summary, sort_keys=True, indent=4, separators=(',', ': '))

#===============================================================================
# MadEvent7CmdShell
#===============================================================================
class MadEvent7CmdShell(MadEvent7Cmd, cmd.CmdShell):
    """The shell command line processor of MadEvent7"""  
    pass
