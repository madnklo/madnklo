#!/usr/bin/env python
from distutils.core import setup, Command
from distutils.command.clean import clean
from distutils.command.install import install
from distutils.command.build_ext import build_ext
from Cython.Build import cythonize

import numpy
import os
pjoin = os.path.join
root_path = os.path.dirname(os.path.realpath( __file__ ))

cython_sources = [
    pjoin(root_path,'madgraph','integrator','phase_space_generators.pyx'),
]

def cythonize_sources():
    return sum([
        cythonize(src_path, compiler_directives={'profile': True},)
        for src_path in cython_sources ],
        [] )

class MG5aMCClean(clean):
    description = "Cleans the build directory and C-sources and shared object libraries coming from Cythonization."

    MG5_aMC_user_options = [
        ('sharedlibraries', None, 'Remove shared object libraries')
    ]
    clean.user_options += MG5_aMC_user_options
    
    def initialize_options(self):
        clean.initialize_options(self)
        self.sharedlibraries = False 
        self.to_clean        = cython_sources

    def finalize_options(self):
        clean.finalize_options(self)
        self.sharedlibraries = bool(self.sharedlibraries)
        assert isinstance(self.to_clean, list) and all(isinstance(p, str) for p in self.to_clean), 'self.to_clean must be a list of paths.'
        assert self.sharedlibraries in [None, False, 1], 'Option shared_libraries must be a boolean, not %s.'%self.sharedlibraries

    def run(self):
        self.all = True
        self.finalize_options()
        clean.run(self)

        # Now clean the C-source and shared object library
        for path in self.to_clean:
            if not path.endswith('.pyx'):
                continue
            c_path = path[:-4]+'.c'
            if os.path.isfile(c_path):
                os.remove(c_path)
            if self.sharedlibraries:
                so_path = c_path = path[:-4]+'.so'
                if os.path.isfile(so_path):
                    os.remove(so_path)


class MG5aMCCompile(Command):
    description = "Compiles MG5aMC Cython sources"
    user_options = [
        ('clean', None, 'Clean up the build directory after compilation'),
        ('force=', None, 'Force recompilation')
    ]

    def initialize_options(self):
        self.clean  = False 
        self.force  = None 

    def finalize_options(self):
        self.clean = bool(self.clean)
        assert isinstance(self.clean, bool), 'Option clean must be a boolean, not %s.'%self.clean
        assert self.force is None or isinstance(self.force, str), "Option force must be 'all' or a specific relative MG5aMC .pyx path."

    def run(self):

        # Force a clean-up if necessary
        if self.force:
            clean_distrib = MG5aMCClean(self.distribution)
            clean_distrib.shared_libraries = True
            if self.force!='all':
                clean_distrib.to_clean = [pjoin(root_path,self.force),]
            clean_distrib.finalize_options()
            clean_distrib.run()
        
        # Assign the cythonized sources (and regenerate the C-sources if not present anymore
        self.distribution.ext_modules = cythonize_sources()
        self.distribution.include_dirs = [numpy.get_include()]

        buildext = build_ext(self.distribution)
        buildext.inplace = 1
        buildext.finalize_options()
        buildext.run()
        if self.clean:
            clean_distrib = MG5aMCClean(self.distribution)
            clean_distrib.shared_libraries = False
            clean_distrib.finalize_options()
            clean_distrib.run()

setup(
    # This attributes will be set directly during the run of MG5aMCCompile
    ext_modules = [],
    cmdclass = {
        'compile' : MG5aMCCompile,
        'clean'   : MG5aMCClean
    }, requires=['scipy', 'matplotlib', 'numpy']
)
