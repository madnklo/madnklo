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
import madgraph.integrator.mappings as mappings
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
            'limit_type'              : None,
            'limit_pattern'           : None,
            'walker'                  : None,
            'process'                 : {'in_pdgs'  : None,
                                       'out_pdgs' : None,
                                       'n_loops'  : None},
            'seed'                    : None,
            'n_steps'                 : 30,
            'min_scaling_variable'    : 1.0e-6,
            'acceptance_threshold'    : 1.0e-6,
            'compute_only_limit_defining_counterterm' : False,
            'include_all_flavors'     : False,
            'apply_higher_multiplicity_cuts' : True,
            'apply_lower_multiplicity_cuts'  : True
        }

        if mode=='poles':
            # Remove some options for the pole
            del testlimits_options['limit_type']
            del testlimits_options['limit_pattern']
            del testlimits_options['walker']
            del testlimits_options['n_steps']
            del testlimits_options['min_scaling_variable']
            del testlimits_options['compute_only_limit_defining_counterterm']
            del testlimits_options['apply_lower_multiplicity_cuts']
        elif mode=='limits':
            del testlimits_options['include_all_flavors']            
 
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
            elif key in ['--compute_only_limit_defining_counterterm','--o'] and mode=='limits':
                if value is None:
                    testlimits_options['compute_only_limit_defining_counterterm'] = True
                else:
                    try:
                        testlimits_options['compute_only_limit_defining_counterterm'] = bool(eval(value))
                    except:
                        raise InvalidCmd("'%s' is not a valid float for option '%s'"%(value, key))                        
            elif key in ['--acceptance_threshold']:
                try:
                    testlimits_options[key[2:]] = float(value)                  
                except ValueError:
                    raise InvalidCmd("'%s' is not a valid float for option '%s'"%(value, key))
            elif key in ['--min_scaling_variable'] and mode=='limits':
                try:
                    testlimits_options[key[2:]] = float(value)                  
                except ValueError:
                    raise InvalidCmd("'%s' is not a valid float for option '%s'"%(value, key))
            elif key in ['--subtraction_mappings_scheme', '--walker'] and mode=='limits':
                try:
                    if value == "None":
                        testlimits_options['walker'] = None
                    else:
                        testlimits_options['walker'] = \
                            mappings.walker_classes_map[value]
                except KeyError:
                    raise InvalidCmd("'%s' is not a valid %s" % (value, key[2:]))
            elif key in ['--include_all_flavors', '--af'] and mode=='poles':
                if value is None:
                    testlimits_options['include_all_flavors'] = True
                else:
                    try:
                        testlimits_options['include_all_flavors'] = \
                                                 {'true':True,'false':False}[value.lower()]
                    except:
                        raise InvalidCmd("'%s' is not a valid float for option '%s'"%(value, key))
            elif key in ['--limit_type','--lt'] and mode=='limits':
                if not isinstance(value, str):
                    raise InvalidCmd("'%s' is not a valid option for '%s'"%(value, key))
                testlimits_options['limit_type'] = value
                if value.lower() == 'soft':
                    testlimits_options['limit_pattern'] = re.compile(r'.*S.*')
                elif value.lower() == 'collinear':
                    testlimits_options['limit_pattern'] = re.compile(r'.*C.*')
                elif value.lower() == 'all':
                    testlimits_options['limit_pattern'] = re.compile(r'.*')
                else:
                    if any(value.startswith(start) for start in ['r"',"r'"]):
                        testlimits_options['limit_pattern'] = re.compile(value)
                    else:
                        # If not specified as a raw string,
                        # we take the liberty of adding the enclosing parenthesis.
                        if not value.startswith('('):
                            value = '(%s,)'%value
                        # If the specified re was not explicitly made a raw string,
                        # we take the liberty of escaping the parenthesis
                        # since this is presumably what the user expects.
                        testlimits_options['limit_pattern'] = re.compile(
                            value.replace('(','\(').replace(')','\)'))
                        
            elif key == '--seed':
                try:
                    testlimits_options['seed'] = int(value)
                except ValueError:
                    raise InvalidCmd("Cannot set '%s' option to '%s'."%(key, value))
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

            else:
                raise InvalidCmd("Option '%s' for the test_IR_%s command not recognized."%(key,mode))        
        
        return new_args, testlimits_options
    
    def parse_launch(self, args):
        """ Parses the argument of the launch command."""

        launch_options = {'integrator': 'VEGAS3',
                          'n_points': None,
                          'n_iterations':None,
                          'verbosity':1,
                          'refresh_filters':'auto',
                          'compile':'auto',
                          'batch_size': 1000,
                          # Here we store a list of lambda function to apply as filters
                          # to the ingegrand we must consider
                          'integrands': [lambda integrand: True]}        
        
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
            elif key in ['--batch_size','--bs']:
                try:
                    launch_options['batch_size'] = int(value)
                except ValueError:
                    raise InvalidCmd("Value '%s' not valid integer for option '%s'."%(value,key))                
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
            
            else:
                raise InvalidCmd("Option '%s' for the launch command not reckognized."%key)

        return launch_options

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
                raise InvalidCmd("Option '%s' for the display integrands command not reckognized."%key)

        return display_options

    def get_integrand_filters(self, filters, mode):
        """ Given user defined filters (a string), returns a lambda function on
        a ME7 integrand that applies it. mode can either be 'accept' or 'reject',
        depending on wheter one wants to keep or reject those integrands."""
        
        assert(mode in ['select','reject'])
        
        if filters is None:
            InvalidCmd("Option '--integrands' must specify values with '='.")
        
        expansion_orders = []
        integrand_types  = []
        
        integrand_type_short_cuts = dict( (k, ( eval('ME7_integrands.ME7Integrand_%s'%k), )
                                            ) for k in ['R','RR','RRR','V','LIB','B'] )
        
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
                    raise InvalidCmd("Option '%s' cannot be "%filter+
                                              "understood as an integrand type specifier.")
        
        filter_functions = []
        if mode == 'select':
            # First the ordering filter if specified
            if expansion_orders:
                filter_functions.append(lambda integrand: 
                    integrand.contribution_definition.correction_order in expansion_orders)
            # The the expansion orders
            if integrand_types:
                filter_functions.append(lambda integrand: 
                                            isinstance(integrand, tuple(integrand_types) ))
        else:
            # First the ordering filter if specified
            if expansion_orders:
                filter_functions.append(lambda integrand: 
                    not (integrand.contribution_definition.correction_order in expansion_orders) )
            # The the expansion orders
            if integrand_types:
                filter_functions.append(lambda integrand: 
                                        not isinstance(integrand, tuple(integrand_types) ))
        
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
            ME7_dump['all_MEAccessors'], root_path=self.me_dir, model=self.model)

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
        
        if MPI_RANK==0:
            # Create a run output directory
            run_output_path = pjoin(self.me_dir,'Results','run_%s'%self.run_card['run_tag'])
            suffix=''
            suffix_number = 0
            while os.path.exists(run_output_path+suffix):
                suffix_number += 1
                suffix = '_%d'%suffix_number
            run_output_path = run_output_path+suffix
            os.makedirs(run_output_path)

        # Setup parallelization
        self.configure_run_mode(self.options['run_mode'])

        # Re-initialize the integrator
        integrator_name = launch_options['integrator']
        integrator_options = self._integrators[integrator_name][1]
        
        integrator_options['verbosity'] = launch_options['verbosity']
        
        if launch_options['n_points']:
            if integrator_name=='VEGAS3':
                integrator_options['survey_n_points'] = launch_options['n_points']
                integrator_options['refine_n_points'] = 5*launch_options['n_points']
                integrator_options['parallelization'] = self.cluster
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

        if integrator_name=='VEGAS3':        
            integrator_options['cluster'] = self.cluster
            integrator_options['batch_size'] = launch_options['batch_size']
        
        integrands_to_consider = ME7_integrands.ME7IntegrandList([ itg for itg in self.all_integrands if
                           all(filter(itg) for filter in launch_options['integrands']) ])

        self.integrator = self._integrators[integrator_name][0](integrands_to_consider, **integrator_options)
        
        if len(set([len(itgd.get_dimensions()) for itgd in integrands_to_consider]))>1 and integrator_name not in ['NAIVE']:
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
        
        # Wrap the call in a propice environment for the run
        with ME7RunEnvironment( silence = False, accessor_optimization = True, loggers = logger_level ):
            xsec, error = self.integrator.integrate()

        logger.info("="*100)
        logger.info('{:^100}'.format("Cross-section with integrator '%s':"%self.integrator.get_name()),'$MG:color:GREEN')
        logger.info('{:^100}'.format("%.5e +/- %.2e [pb]"%(xsec, error)),'$MG:color:BLUE')
        logger.info("="*100+"\n")
        
        if MPI_RANK==0:
            # Write the result in 'cross_sections.dat' of the result directory
            xsec_summary = open(pjoin(run_output_path,'cross_sections.dat'),'w')
            xsec_summary_lines = []        
            xsec_summary_lines.append('%-30s%-30s%-30s'%('','Cross-section [pb]','MC uncertainty [pb]'))
            xsec_summary_lines.append('%-30s%-30s%-30s'%('Total','%.8e'%xsec,'%.3e'%error))
            xsec_summary.write('\n'.join(xsec_summary_lines))
            xsec_summary.close()
    
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

        for integrand in self.all_integrands:
            if not hasattr(integrand, 'test_IR_limits'):
                continue
            logger.debug(
                'Now testing IR limits of the following integrand:\n' +
                integrand.nice_string() )
            integrand.test_IR_limits(test_options=testlimits_options)

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

        for integrand in self.all_integrands:
            if not hasattr(integrand, 'test_IR_poles'):
                continue
            logger.debug('Now testing IR poles of the following integrand:\n%s'%(integrand.nice_string()))
            integrand.test_IR_poles(test_options=testpoles_options)

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
