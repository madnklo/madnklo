##########################################################################################
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
##########################################################################################

"""Tools to facilitate debugging"""

import ipdb
import madgraph.various.misc as misc

#=====================================
# Function decorators
#=====================================

def watch_errors(f):
    """Start ipdb in the function f if an exception is raised"""
    def wrapped_f(*args,**kwargs):
        with ipdb.launch_ipdb_on_exception():
            return f(*args,**kwargs)
    return wrapped_f

def signal_before(f):
    """Print signal 'PING' before function execution"""
    def wrapped_f(*args,**kwargs):
        misc.sprint("PING")
        return f(*args,**kwargs)
    return wrapped_f

def signal_after(f):
    """Print signal 'PONG' after function execution"""
    def wrapped_f(*args,**kwargs):
        result = f(*args,**kwargs)
        misc.sprint("PING")
        return result
    return wrapped_f

def print_before(msg="PING"):
    """Print message before function execution (default: PING)"""
    def wrapper(f):
        def wrapped_f(*args,**kwargs):
            misc.sprint(msg)
            return f(*args,**kwargs)
        return wrapped_f
    return wrapper

def print_after(msg="PONG"):
    """Print message after function execution (default: PONG)"""
    def wrapper(f):
        def wrapped_f(*args,**kwargs):
            result = f(*args, **kwargs)
            misc.sprint(msg)
            return result
        return wrapped_f
    return wrapper


#=====================================
# Debugger interaction
#=====================================

set_trace = ipdb.set_trace
