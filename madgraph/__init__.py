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
from __future__ import absolute_import
class MadGraph5Error(Exception):
    """Exception raised if an exception is find 
    Those Types of error will stop nicely in the cmd interface"""

class InvalidCmd(MadGraph5Error):
    """a class for the invalid syntax call"""

class aMCatNLOError(MadGraph5Error):
    """A MC@NLO error"""

import os
import logging
import time

#################
# numpy imports
#
import numpy as np
# "cimport" is used to import special compile-time information
# about the numpy module (this is stored in a file numpy.pxd which is
# currently part of the Cython distribution).
##cimport numpy as np
# We now need to fix a datatype for our arrays. I've used the variable
# DTYPE for this, which is assigned to the usual NumPy runtime
# type info object.
DTYPE = np.float64
# "ctypedef" assigns a corresponding compile-time type to DTYPE_t. For
# every type in the numpy module there's a corresponding compile-time
# type with a _t-suffix.
##ctypedef np.float64_t DTYPE_t
#
#################

#Look for basic file position MG5DIR and MG4DIR
MG5DIR = os.path.realpath(os.path.join(os.path.dirname(__file__),
                                                                os.path.pardir))
if ' ' in MG5DIR:
   logging.critical('''\033[1;31mpath to MG5: "%s" contains space. 
    This is likely to create code unstability. 
    Please consider changing the path location of the code\033[0m''' % MG5DIR)
   time.sleep(1)
MG4DIR = MG5DIR
ReadWrite = os.access(MG5DIR, os.W_OK) # W_OK is for writing

if ReadWrite:
    # Temporary fix for problem with auto-update
    try:
        tmp_path = pjoin(MG5DIR, 'Template','LO','Source','make_opts')
        #1480375724 is 29/11/16
        if os.path.exists(tmp_path) and os.path.getmtime(tmp_path) < 1480375724:
            os.remove(tmp_path)
            shutil.copy(pjoin(MG5DIR, 'Template','LO','Source','.make_opts'),
                    pjoin(MG5DIR, 'Template','LO','Source','make_opts'))
    except Exception as error:
        pass

# Check MPI configuration
MPI_ACTIVE = False
MPI_RANK   = 0
MPI_SIZE   = 1
try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    MPI_RANK = comm.Get_rank()
    MPI_SIZE = comm.Get_size()
    if MPI_SIZE > 1:
       MPI_ACTIVE = True
    else:
       MPI_ACTIVE = False
except:
    pass

  
ADMIN_DEBUG = False  
if os.path.exists(os.path.join(MG5DIR,'bin', 'create_release.py')):
    if os.path.exists(os.path.join(MG5DIR,'.bzr')):
        ADMIN_DEBUG = True

if __debug__ or ADMIN_DEBUG:
    ordering = True
else:
    ordering = False
        
