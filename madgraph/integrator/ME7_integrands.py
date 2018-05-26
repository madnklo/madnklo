##########################################################################################
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
##########################################################################################
"""Collection of ME7 Integrand classes."""

import numpy as np
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

try:
    import pyjet
    PYJET_AVAILABLE = True
except:
    PYJET_AVAILABLE = False    

# useful shortcut
pjoin = os.path.join

#root_path = pjoin(os.path.dirname(os.path.realpath( __file__ )), *([os.path.pardir]*3))
#sys.path.insert(0, root_path)

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
import madgraph.integrator.mappings as mappings
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
                     topologies_to_processes,
                     processes_to_topologies,
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
            elif contribution_definition.correction_order == 'NNNLO':
                if contribution_definition.n_loops == 0 and \
                   contribution_definition.n_unresolved_particles == 3:
                    target_type = 'TripleReals'            
                else:
                    raise MadGraph5Error("Some NNNLO type of integrands are not implemented yet.") 
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
                       topologies_to_processes,
                       processes_to_topologies,
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
                                       'topologies_to_processes'    : topologies_to_processes,
                                       'processes_to_topologies'    : processes_to_topologies,
                                       'all_MEAccessors'            : None,
                                       'ME7_configuration'          : None,
                                       'options'                    : {} }
        
        # The original ContributionDefinition instance at the origin this integrand 
        self.contribution_definition    = contribution_definition
        # The process map of the Contribution instance at the origin of this integrand.
        # The format is identical to the one generated from the function 'get_process_map' of a contribution.
        self.processes_map              = processes_map
        
        # Add information about the topology of the diagrams constituting the processes,
        # so as to be able to build efficient phase-space parametrizations. The format of these dictionaries
        # is specified in the function 'set_phase_space_topologies' of the class contributions.Contribution
        self.topologies_to_processes    = topologies_to_processes
        self.processes_to_topologies    = processes_to_topologies
        
        # An instance of accessors.MEAccessorDict providing access to all ME available as part of this
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

    def get_additional_nice_string_printout_lines(self):
        """ Additional printout information for nice_string. 
        Meant to possibly be overloaded by daughter classes."""
        return []

    def get_nice_string_process_line(self, process_key, defining_process, format=0):
        """ Return a nicely formated process line for the function nice_string of this 
        contribution. Can be overloaded by daughter classes."""
        GREEN = '\033[92m'
        ENDC = '\033[0m'        
        return GREEN+'  %s'%defining_process.nice_string(print_weighted=False).\
                                                               replace('Process: ','')+ENDC

    def nice_string(self, format=0):
        """ Nice string representation of self."""
        BLUE = '\033[94m'
        GREEN = '\033[92m'
        ENDC = '\033[0m'
        res = ['%-30s:   %s'%('ME7Integrand_type',type(self))]
        res.extend([self.contribution_definition.nice_string()])
        if not self.topologies_to_processes is None:
            res.append('%-30s:   %d'%('Number of topologies', 
                                                    len(self.topologies_to_processes.keys())))
        res.extend(self.get_additional_nice_string_printout_lines())

        if format < 1:
            res.append('Generated and mapped processes for this contribution: %d (+%d mapped)'%
                       ( len(self.processes_map.keys()),
                         len(sum([v[1] for v in self.processes_map.values()],[])) ) )
        else:
            res.append('Generated and mapped processes for this contribution:')
            for process_key, (defining_process, mapped_processes) in self.processes_map.items():
                res.append(self.get_nice_string_process_line(process_key, defining_process, format=format))                
                for mapped_process in mapped_processes:
                    res.append(BLUE+u'   \u21b3  '+mapped_process.nice_string(print_weighted=False)\
                                                                        .replace('Process: ','')+ENDC)
            
        return '\n'.join(res).encode('utf-8')

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
            # Ask for logging.DEBUG-1 so as to only have lhapdf verbose if really desired.
            lhapdf.setVerbosity(1 if logger.level<=(logging.DEBUG-1) else 0)

            pdfsets_dir = subprocess.Popen([lhapdf_config,'--datadir'],\
                                           stdout=subprocess.PIPE).stdout.read().strip()
            lhapdf.pathsPrepend(pdfsets_dir)
            lhapdf_version = subprocess.Popen([lhapdf_config,'--version'],\
                                           stdout=subprocess.PIPE).stdout.read().strip()
            pdf_info = common_run.CommonRunCmd.get_lhapdf_pdfsets_list_static(pdfsets_dir, lhapdf_version)
            lhaid = self.run_card.get_lhapdf_id()
            if lhaid not in pdf_info:
                raise InvalidCmd("Could not find PDF set with lhaid #%d in %s."%(lhaid, pdfsets_dir))
            pdf_set_name = pdf_info[lhaid]['filename']
            if not os.path.isdir(pjoin(pdfsets_dir, pdf_set_name)):
                raise InvalidCmd("Could not find PDF set directory named "+
                    "'%s' in '%s'.\n"%(pdf_set_name, pdfsets_dir)+
                    "It can be downloaded from LHAPDF official online resources.")

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
        
        return cls(  model, 
                     run_card,                             
                     dump['contribution_definition'],
                     dump['processes_map'],
                     dump['topologies_to_processes'],
                     dump['processes_to_topologies'],
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

    def is_part_of_process_selection(self, process_list, selection=None):
        """ Checks whether any of the specified processes in the process_list provided matches the user's process
        selection. If not provided, returns True by default. 'selection' is a dictionary with the format:
           {'in_pdgs'  : ( (in_pgs1), (in_pdgs2), ...)
            'out_pdgs' : ( (out_pgs1), (out_pdgs2), ...) 
            'n_loops'  : n_loops }"""
        
        def pdg_list_match(target_list, selection_list):
            if len(target_list) != len(selection_list):
                return False
            targets = dict( (k, target_list.count(k)) for k in set(target_list) )
            found_it = False
            for sel in itertools.product(*selection_list):
                found_it = True
                for k, v in targets.items():
                    if sel.count(k) != v:
                        found_it = False
                        break
                if found_it:
                    break
            return found_it                
    
        for process in process_list:

            if (not selection['in_pdgs'] is None) and \
               (not pdg_list_match(process.get_initial_ids(), selection['in_pdgs'])):
                continue

            if (not selection['out_pdgs'] is None) and \
               (not pdg_list_match(process.get_final_ids_after_decay(), selection['out_pdgs'])):
                continue
            
            if (not selection['n_loops'] is None) and process.get('n_loops') != selection['n_loops']:
                continue
            return True

        return False

    def find_counterterms_matching_regexp(
        self, counterterms, limit_pattern=None):
        """ Find all mappings that match a particular limits given in argument
        (takes a random one if left to None). This function is placed here given that
        it can be useful for both the ME7Integrnd_V and ME7_integrand_R."""

        # First select only the counterterms which are not pure matrix elements 
        # (i.e. they have singular structures) and also exclude here soft-collinear 
        # counterterms since they do not have an approach limit function.
        selected_counterterms = [ ct for ct in counterterms if ct.is_singular() ]
        if len(selected_counterterms)==0:
            return []
        
        limit_pattern_re = None
        if limit_pattern:
            if limit_pattern.lower() == 'soft':
                limit_pattern_re = re.compile(r'.*S.*')
            elif limit_pattern.lower() == 'collinear':
                limit_pattern_re = re.compile(r'.*C.*')
            elif limit_pattern.lower() == 'puresoft':
                limit_pattern_re = re.compile(r'^[S\d,\(\)]*$')
            elif limit_pattern.lower() == 'purecollinear':
                limit_pattern_re = re.compile(r'^[C\d,\(\)]*$')
            elif limit_pattern.lower() == 'all':
                limit_pattern_re = re.compile(r'.*')
            elif any(limit_pattern.startswith(start) for start in ['r"', "r'"]):
                limit_pattern_re = re.compile(eval(limit_pattern))
            else:
                # Check if a list of counterterms is specified
                try:
                    list_limit_pattern = eval(limit_pattern)
                    if not isinstance(list_limit_pattern, list):
                        raise BaseException
                except:
                    list_limit_pattern = [limit_pattern]
                new_list_limit_pattern = []
                for limit_pattern in list_limit_pattern:
                    if not limit_pattern.startswith('('):
                        # If not specified as a raw string, we take the liberty of adding 
                        # the enclosing parenthesis.
                        limit_pattern = '(%s,)' % limit_pattern
                    # We also take the liberty of escaping the parenthesis
                    # since this is presumably what the user expects.
                    limit_pattern = limit_pattern.replace('(', '\(').replace(')', '\)')
                    new_list_limit_pattern.append(limit_pattern)
                limit_pattern_re = re.compile(r'^(%s)$'%(
                    '|'.join(limit_pattern for limit_pattern in new_list_limit_pattern) ))    

        returned_counterterms = []
        if not limit_pattern:
            returned_counterterms.append(random.choice(selected_counterterms))
        else:
            for counterterm in selected_counterterms:
                if re.match(limit_pattern_re, str(counterterm)):
                    returned_counterterms.append(counterterm)

        return returned_counterterms

    def pass_flavor_blind_cuts(self, PS_point, process_pdgs, n_jets_allowed_to_be_clustered = None):
        """ Implementation of a minimal set of isolation cuts. This can be made much nicer in the future and 
        will probably be taken outside of this class so as to ease user-specific cuts, fastjet, etc...
        We consider here a two-level cuts system, this first one of which is flavour blind.
        The 'n_jets_allowed_to_be_clustered' is an option that allows to overwrite the
        maximum number of jets that can be clustered and which is by default taken to be:
            self.contribution_definition.n_unresolved_particles 
        This is useful when using this function to apply cuts to the reduced PS of the CTs."""
        debug_cuts = False
        # This is a temporary function anyway which should eventually be replaced by a full
        # fledged module for handling generation level cuts, which would also make use of fjcore.
        # return True

        from madgraph.integrator.vectors import LorentzVectorDict, LorentzVectorList, LorentzVector
        
        # These cuts are not allowed to resolve flavour, but only whether a particle is a jet or not
        def is_a_jet(pdg):
            return abs(pdg) in range(1,self.run_card['maxjetflavor']+1)+[21]
        
        def is_a_lepton(pdg):
            return abs(pdg) in [11,13,15]
        
        def is_a_neutrino(pdg):
            return abs(pdg) in [12,14,16]
        
        def is_a_photon(pdg):
            return pdg==22
        
        if debug_cuts: logger.debug( "Processing flavor-blind cuts for process %s and PS point:\n%s"%(
            str(process_pdgs), LorentzVectorList(PS_point).__str__(n_initial=self.phase_space_generator.n_initial) ))

        if n_jets_allowed_to_be_clustered is None:
            n_jets_allowed_to_be_clustered = self.contribution_definition.n_unresolved_particles

        ###################################################################################
        # JET CLUSTERING AND CUTS
        ###################################################################################
        
        ptj_cut = self.run_card['ptj']
        drjj_cut = self.run_card['drjj']
        etaj_cut = self.run_card['etaj']
        
        if ptj_cut <= 0. and    \
           drjj_cut <= 0. and   \
           etaj_cut <= 0. :
            return True
        else:
            # If fastjet is needed but not found, make sure to stop
            if (not PYJET_AVAILABLE) and n_jets_allowed_to_be_clustered>0:
                raise MadEvent7Error("Fast-jet python bindings are necessary for integrating"+
                             " real-emission type of contributions. Please install pyjet.")

        # The PS point in input is sometimes provided as a dictionary or a flat list, but
        # we need it as a flat list here, so we force the conversion
        if isinstance(PS_point, LorentzVectorDict):
            PS_point = PS_point.to_list()
        elif not isinstance(PS_point, LorentzVectorList):
            PS_point = LorentzVectorList(LorentzVector(v) for v in PS_point)    

        if PYJET_AVAILABLE and drjj_cut > 0.:

            # First identify all partonic jets
            jets_list = []
            for i, p in enumerate(PS_point[self.n_initial:]):
                if is_a_jet(process_pdgs[1][i]):
                    jets_list.append(tuple(list(p)+[i+self.n_initial+1,]))
            # Count partonic jets
            starting_n_jets = len(jets_list)

            # Cluster them with fastjet
            this_event = np.array(jets_list,dtype=np.dtype(
                    [('E', 'f8'), ('px', 'f8'), ('py', 'f8'), ('pz', 'f8'), ('id', 'i8')]) )
            sequence = pyjet.cluster(this_event, R=drjj_cut, p=-1, ep=True)
            jets = sequence.inclusive_jets()

#            misc.sprint('Process flavors: %s'%str(process_pdgs))
#            misc.sprint('Starting from configuration:\n%s'%str(PS_point))
#            misc.sprint('The following %d jets were found:'%len(jets))
#            misc.sprint("{0: <5} {1: >10} {2: >10} {3: >10} {4: >10} {5: >10} {6: >10} {7: >10}".format(
#                                      "jet#", "E", "Px", "Py", "P.z", "mass", "pT", "#constit."))
#            for i, jet in enumerate(jets):
#                misc.sprint("{0: <5} {1: >10} {2: >10} {3: >10} {4: >10} {5: >10} {6: >10} {7: >10}".format(
#                    i + 1, jet.e, jet.px, jet.py, jet.pz, jet.mass, jet.pt, 
#                    str(tuple(constituent.id for constituent in jet))
#                    ))
#                misc.sprint('\n'+'\n'.join('  | %s'%constituent for
#                                           constituent in jet.constituents_array(ep=True)))

            # Remove jets whose pT is below the user defined cut:
            if ptj_cut > 0.:
                jets = [jet for jet in jets if jet.pt >= ptj_cut]
            
#            # A useless (in production) consistency check for deltaR:
#            for i, j1 in enumerate(jets):
#                for j, j2 in enumerate(jets):
#                    if i==j: continue
#                    pj1 = LorentzVector([j1.e, j1.px, j1.py, j1.pz])
#                    pj2 = LorentzVector([j2.e, j2.px, j2.py, j2.pz])
#                    if pj1.deltaR(pj2) < drjj_cut:
#                        raise MadGraph5Error("Inconsistency with fastjet in pass_cuts : "+
#                                                  " %.5f < %.5f"%(pj1.deltaR(pj2), drjj_cut))
#                    misc.sprint(pj1.deltaR(pj2), drjj_cut)                
            
            # Make sure that the number of clustered jets is at least larger or equal to the
            # starting list of jets minus the number of particles that are allowed to go
            # unresolved in this contribution.
            if debug_cuts: logger.debug("Number of identified jets: %d (min %d)"%
                           ( len(jets), (starting_n_jets-n_jets_allowed_to_be_clustered) ))
            if len(jets) < (starting_n_jets-n_jets_allowed_to_be_clustered):

                return False
            
            all_jets = LorentzVectorList([LorentzVector(
                               [a_jet.e, a_jet.px, a_jet.py, a_jet.pz]) for a_jet in jets])
            
        else:
            all_jets = LorentzVectorList([p for i, p in enumerate(PS_point[self.n_initial:])
                                                          if is_a_jet(process_pdgs[1][i])])
            if ptj_cut > 0.:
                # Apply the Ptj cut first
                for i, p in enumerate(all_jets):
                    if debug_cuts: logger.debug('pj_%i.pt()=%.5e'%((i+1),p.pt()))
                    if p.pt() < ptj_cut:
                        return False
    
            if drjj_cut > 0.:
                # And then the drjj cut
                for i, p1 in enumerate(all_jets):
                    for j, p2 in enumerate(all_jets):
                        if j <= i:
                            continue
                        if debug_cuts: logger.debug('deltaR(pj_%i,pj_%i)=%.5e'%(
                             i+1, j+1, p1.deltaR(p2)))
                        if p1.deltaR(p2) < drjj_cut:
                            return False
    
        # Now handle all other cuts
        if etaj_cut > 0.:
            for i, p_jet in enumerate(all_jets):
                if debug_cuts: logger.debug('eta(pj_%i)=%.5e'%(i+1,p_jet.pseudoRap()))
                if abs(p_jet.pseudoRap()) > etaj_cut:
                    return False

        ###################################################################################
        # LEPTON AND PHOTON CUTS
        ###################################################################################
        
        for i, p in enumerate(PS_point[self.n_initial:]):
            # photons
            if is_a_photon(process_pdgs[1][i]):
                if debug_cuts: logger.debug('pta_%i.pt()=%.5e'%((i+1),p.pt()))
                if self.run_card['pta'] > 0.0 and p.pt() < self.run_card['pta']:
                    return False
                if debug_cuts: logger.debug('eta(pa_%i)=%.5e'%(i+1,p.pseudoRap()))
                if self.run_card['etaa'] > 0.0 and abs(p.pseudoRap()) > self.run_card['etaa']:
                    return False
                for j, p2 in enumerate(PS_point[self.n_initial:]):
                    if j > i and is_a_photon(process_pdgs[1][j]):
                        if debug_cuts: logger.debug('deltaR(pa_%i,pa_%i)=%.5e'%(i+1, j+1, p.deltaR(p2)))
                        if self.run_card['draa'] > 0.0 and p.deltaR(p2) < self.run_card['draa']:
                            return False
                    if is_a_lepton(process_pdgs[1][j]):
                        if debug_cuts: logger.debug('deltaR(pa_%i,pl_%i)=%.5e'%(i+1, j+1, p.deltaR(p2)))
                        if self.run_card['dral'] > 0.0 and p.deltaR(p2) < self.run_card['dral']:
                            return False
                for j, p_jet in enumerate(all_jets):
                    if debug_cuts: logger.debug('deltaR(pa_%i,pj_%i)=%.5e'%(i+1, j+1, p.deltaR(p_jet)))
                    if self.run_card['draj'] > 0.0 and p.deltaR(p_jet) < self.run_card['draj']:
                        return False
            
            # leptons
            if is_a_lepton(process_pdgs[1][i]):
                if debug_cuts: logger.debug('ptl_%i.pt()=%.5e'%((i+1),p.pt()))
                if self.run_card['ptl'] > 0.0 and p.pt() < self.run_card['ptl']:
                    return False
                if debug_cuts: logger.debug('eta(pl_%i)=%.5e'%(i+1,p.pseudoRap()))
                if self.run_card['etal'] > 0.0 and abs(p.pseudoRap()) > self.run_card['etal']:
                    return False
                for j, p2 in enumerate(PS_point[self.n_initial:]):
                    if j <= i:
                        continue
                    if is_a_lepton(process_pdgs[1][j]):
                        if debug_cuts: logger.debug('deltaR(pl_%i,pl_%i)=%.5e'%(i+1, j+1, p.deltaR(p2)))
                        if self.run_card['drll'] > 0.0 and p.deltaR(p2) < self.run_card['drll']:
                            return False
                for j, p_jet in enumerate(all_jets):
                    if debug_cuts: logger.debug('deltaR(pl_%i,pj_%i)=%.5e'%(i+1, j+1, p.deltaR(p_jet)))
                    if self.run_card['drjl'] > 0.0 and p.deltaR(p_jet) < self.run_card['drjl']:
                        return False  
            

        # All cuts pass, therefore return True
        return True

    def pass_flavor_sensitive_cuts(self, PS_point, flavors):
        """ Implementation of a minimal set of isolation cuts. This can be made much nicer in the future and 
        will probably be taken outside of this class so as to ease user-specific cuts, fastjet, etc...
        We consider here a two-level cuts system, this second one of which is flavour sensitive."""

        # None implemented yet
        return True
    
        debug_cuts = False
        from madgraph.integrator.vectors import LorentzVectorDict, LorentzVectorList, LorentzVector

        # The PS point in input is sometimes provided as a dictionary or a flat list, but
        # we need it as a flat list here, so we force the conversion
        if isinstance(PS_point, LorentzVectorDict):
            PS_point = PS_point.to_list()
        elif not isinstance(PS_point, LorentzVectorList):
            PS_point = LorentzVectorList(LorentzVector(v) for v in PS_point)   

        if debug_cuts: logger.debug( "Processing flavor-sensitive cuts for flavors %s and PS point:\n%s"%(
            str(flavors), LorentzVectorList(PS_point).__str__(n_initial=self.phase_space_generator.n_initial) ))
        
        # All cuts pass, therefore return True
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
        # /!\ WARNING One typically only boosts to the lab frame at the very end. But we do it here nonetheless. Can easily be changed.
        ###
        # self.phase_space_generator.boost_to_lab_frame(PS_point, xb_1, xb_2)
                
        if __debug__: logger.debug("Considering the following PS point:\n%s"%(PS_point.__str__(
                                            n_initial=self.phase_space_generator.n_initial) ))
        
        # Account for PS weight
        wgt *= PS_weight
        if __debug__: logger.debug("PS_weight: %.5e"%PS_weight)
        
        # Include the flux factor
        flux = 1.
        if self.n_initial == 2:
            flux = 1. / (2.*math.sqrt(self.Lambda(E_cm**2, self.masses[0][0]**2, self.masses[0][1]**2)))
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
        for process_key, (process, mapped_processes) in self.processes_map.items():
            if __debug__: logger.debug('Now considering the process group from %s.'%process.nice_string())
            
            this_process_wgt = wgt
            process_pdgs = process.get_cached_initial_final_pdgs()

            # Apply flavor blind cuts
            if not self.pass_flavor_blind_cuts(PS_point, process_pdgs):
                if __debug__: logger.debug('Event failed the flavour_blind generation-level cuts.')
                if __debug__: logger.debug(misc.bcolors.GREEN + 'Returning a weight of 0. for this integrand evaluation.' + misc.bcolors.ENDC)            
                return 0.0
            
            all_processes = [process,]+mapped_processes
            all_process_pdgs = []
            # Add mirror processes if present
            for proc in all_processes:
                initial_final_pdgs = proc.get_cached_initial_final_pdgs()
                all_process_pdgs.append(initial_final_pdgs)
                if proc.get('has_mirror_process'):
                    all_process_pdgs.append(
                       ( tuple(list(initial_final_pdgs[0])[::-1]), initial_final_pdgs[1] ) )

            if __debug__:
                if abs(self.run_card['lpp1'])==abs(self.run_card['lpp2'])==1:
                    logger.debug("Bjorken x's x1, x2, sqrt(x1*x2*s): %.5e, %.5e. %.5e"%(
                                 xb_1, xb_2, math.sqrt(self.collider_energy**2*xb_1*xb_2)))
            proc_PDFs_weights = []
            for proc_pdgs in all_process_pdgs:            
                if self.run_card['lpp1']==self.run_card['lpp2']==1:            
                    PDF1 = self.get_pdfQ2(self.pdf, proc_pdgs[0][0], xb_1, mu_f1**2)
                    PDF2 = self.get_pdfQ2(self.pdf, proc_pdgs[0][1], xb_2, mu_f2**2)
                    proc_PDFs_weights.append(PDF1*PDF2)
                    if __debug__: logger.debug("PDF(x1, %d), PDF(x2, %d) = %.5e, %.5e"%(
                                             proc_pdgs[0][0], proc_pdgs[0][1], PDF1, PDF2))
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
            this_process_wgt *= sum(proc_PDFs_weights)
    
            # Finally include the short-distance weight, with the PS point in a dictionary format.
            sigma_wgt = self.sigma(
                PS_point.to_dict(), process_key, process, selected_flavors, this_process_wgt, 
                                                                        mu_r, mu_f1, mu_f2)
            if __debug__: logger.debug('Short-distance sigma weight for this subprocess: %.5e'%sigma_wgt)        
            this_process_wgt *= sigma_wgt
            
            # Accumulate this process weight
            total_wgt += this_process_wgt

        # Now finally return the total weight for this contribution
        if __debug__: logger.debug(misc.bcolors.GREEN + "Final weight returned: %.5e"%total_wgt + misc.bcolors.ENDC)
        if __debug__: logger.debug("="*80)
        
        return total_wgt
    
    def sigma(self, PS_point, process_key, process, flavors, process_wgt, mu_r, mu_f1, mu_f2):
        """ 
        This is the core function of the integrand where the short-distance objects like the matrix elements,
        the counterterms, the mappings, etc.. will be evaluated.
        The process has a lot of extra meta information that sigma may consider, but flavors provides what flavour 
        assignment should be used for some of the flavour dependent steps (like calling the cuts for instance).
        Notice that this function may be called several times with the exact same arguments but different 
        processes / flavour, so some caching may be in order. 
        The output of sigma, sigma_wgt, should be the weight summed over the various contributions building sigma
        (ME's and counterterms).
        Finally, this function is also responsible for calling the observable, using process_wgt*sigma_wgt.
        It can potentially call it several times, for the various counterterms.
        """

        sigma_wgt = 1.

        alpha_s = self.model.get('parameter_dict')['aS']

        # Access to the matrix element. ME_result is an instance of a subclassed dictionary which includes
        # all independent results available as of now for the particular arguments specified.
        # We specify pdgs to None her to avoid the need of any permutation since we follow the order of
        # the defining process here which is the one that was exported.
        # For the reduced matrix elements however, this cannot be done.
        ME_evaluation, all_results = self.all_MEAccessors(
                                            process, PS_point, alpha_s, mu_r, pdgs=flavors)
        
        ## Some code to test the color correlated MEs
        ##color_correlation_to_consider = ( ((3, -2, 3), (-2, -2, -1)), ((3, -1, 3), (-1, -1, -2)) )
        ##color_correlation_to_consider = ( -1, ((3, -1, 3), (-1, -1, -2)) )
        ##color_correlation_to_consider = ( -1, ((3, -1, 3), ) )        
        ##ME_evaluation, all_results = self.all_MEAccessors(
        ##    process, PS_point, alpha_s, mu_r, pdgs=flavors, 
        ##    color_correlation=( color_correlation_to_consider, ),
        ##    return_all_res = True)
        ##misc.sprint(str(ME_evaluation))
        ##misc.sprint(str(all_results))
        ##ME_evaluation, all_results = self.all_MEAccessors(
        ##    process, PS_point, alpha_s, mu_r, pdgs=flavors, 
        ##    color_correlation=( ( ((4, -1, 4), ), ((3, -1, 3), ) ), ),
        ##    return_all_res = True)
        ##misc.sprint(str(ME_evaluation))
        ##misc.sprint(str(all_results))
        
        sigma_wgt *= ME_evaluation['finite']
        
        if self.apply_observables:
            data_for_observables = {'PS_point': PS_point, 'flavors' : flavors}
            self.observable_list.apply_observables(sigma_wgt*process_wgt, data_for_observables)

        return sigma_wgt

##########################################################################################
# Daughter classes for the various type of ME7Integrands which require a redefinition 
# of sigma() and possibly other functions.
##########################################################################################
class ME7Integrand_B(ME7Integrand):
    """ME7Integrand for the computation of a Born type of contribution."""
    
    def sigma(self, PS_point, process_key, process, flavors, process_wgt,
                                                        mu_r, mu_f1, mu_f2, *args, **opts):

        return super(ME7Integrand_B, self).sigma(
            PS_point, process_key, process, flavors, process_wgt,
            mu_r, mu_f1, mu_f2, *args, **opts
        )

class ME7Integrand_LIB(ME7Integrand):
    """ ME7Integrand for the computation of a Loop-Induced Born type of contribution."""
    def sigma(self, PS_point, process_key, process, flavors, process_wgt, mu_r, mu_f1, mu_f2, *args, **opts):
        return super(ME7Integrand_LIB, self).sigma(PS_point, process_key, process, flavors, process_wgt, 
                                                                        mu_r, mu_f1, mu_f2, *args, **opts)
        
class ME7Integrand_V(ME7Integrand):
    """ ME7Integrand for the computation of a one-loop virtual type of contribution."""
    
    def __init__(self, *args, **opts):
        """Initialize a virtual type of integrand, adding additional relevant attributes."""
        
        try:  
            self.integrated_counterterms = opts.pop('integrated_counterterms')
        except KeyError:
            raise MadEvent7Error(
                "Constructor of class ME7Integrand_V requires the option "
                "'integrated_counterterms' to be specified."
            )

        super(ME7Integrand_V, self).__init__(*args, **opts)
        self.initialization_inputs['options']['integrated_counterterms'] = \
                                                              self.integrated_counterterms

    def get_additional_nice_string_printout_lines(self, format=0):
        """ Return additional information lines for the function nice_string of this contribution."""
        res = []
        if self.integrated_counterterms:
            res.append('%-30s:   %d'%('Nb. of integrated counterterms', 
                                      len(sum(self.integrated_counterterms.values(),[]))))
        return res
    
    def get_nice_string_process_line(self, process_key, defining_process, format=0):
        """Return a nicely formatted process line for the function nice_string of this
        contribution."""

        GREEN = '\033[92m'
        ENDC = '\033[0m'        
        res = GREEN+'  %s'%defining_process.nice_string(print_weighted=False).\
                                                               replace('Process: ','')+ENDC

        if not self.integrated_counterterms:
            return res

        if format<2:
            if process_key in self.integrated_counterterms:
                res += ' | %d integrated counterterms'%len(self.integrated_counterterms[process_key])
            else:
                res += ' | 0 integrated counterterm'
                
        else:
            long_res = [' | with the following integrated counterterms:']
            for CT_properties in self.integrated_counterterms[process_key]:
                CT = CT_properties['integrated_counterterm']
                if format==2:
                    long_res.append( '   | %s'%CT.__str__(
                                        print_n=True, print_pdg=False, print_state=False )  )
                elif format==3:
                    long_res.append( '   | %s'%CT.__str__(
                                        print_n=True, print_pdg=True, print_state=True )  )
                elif format==4:
                    long_res.append( '   | %s'%str(CT))
                elif format>4:
                    long_res.append( '   | %s'%str(CT))
                    for key, value in CT_properties.items():
                        if not key in ['integrated_counterterm', 'matching_process_key']:
                            long_res.append( '     + %s : %s'%(key, str(value)))

            res += '\n'.join(long_res)

        return res

    def evaluate_integrated_counterterm(self, integrated_CT_characteristics, 
                    PS_point, flavors, input_mapping, hel_config=None, compute_poles=True):
        """ Evaluates the specified integrated counterterm specified along with its other
        characteristics, like for example the list of flavors assignments that the resolved
        process it corresponds to can take. This function returns
            integrated_CT_res, reduced_flavors
        where integrated_CT_res is the result of the evaluation (and EpsilonExpansion instance
        if compute_poles is True, otherwise just a float, and reduced_flavors is the list
        of flavors for the reduced process to be used when calling the pass_cuts and 
        observable functions.
        """

        # Access the various characteristics of the integrated counterterm passed to this
        # function.
        counterterm = integrated_CT_characteristics['integrated_counterterm'] 
        resolved_flavors = integrated_CT_characteristics['resolved_flavors_combinations']
        reduced_flavors = integrated_CT_characteristics['reduced_flavors_combinations']
        symmetry_factor = integrated_CT_characteristics['symmetry_factor']

        n_initial = len(flavors[0])
        n_final   = len(flavors[1])
        
        mapped_flavors = ( 
            tuple( flavors[0][input_mapping[i]] for i in range(n_initial) ),
            tuple( flavors[1][input_mapping[i]-n_initial] for i in 
                                                      range(n_initial, n_initial+n_final) )
        )        
        # And the multiplicity prefactor coming from the several *resolved* flavor assignment
        # that this counterterm can lead to. Typically an integrated counterterm for g > q qbar
        # splitting will have the same reduced flavors, but with the gluon coming from
        # n_f different massless flavors. So that its multiplcitiy factor is n_f.        
        # Also, if you think of the resolved flavors e+ e- > c c~ u u~, the two counterterms
        # C(3,4) and C(5,6) will lead to the reduced flavors g u u~ and c c~ g respectively.
        # These two are however mapped in the same virtual group, so it's possible that if
        # the flavors 'g u u~' was selected for example, then the integrated CT C(5,6) would
        # not have any match in its reduced flavors which contains 'c c~ g' only, and it should
        # therefore be skipped. 
        if not (mapped_flavors in reduced_flavors):
            return None, None
        integrated_CT_multiplicity = reduced_flavors[mapped_flavors]
        # Also account for the multiplicity from overall symmetry factors S_t
        integrated_CT_multiplicity *= symmetry_factor

        if isinstance(PS_point,dict):
            # Dictionary format LorentzVectorDict starts at 1
            mapped_PS_point = phase_space_generators.LorentzVectorDict(   
                (i+1, PS_point[input_mapping[i]+1]) for i in range(n_initial+n_final) )
        else:
            # List formatLorentzVectorList starts at 0
            mapped_PS_point = phase_space_generators.LorentzVectorDict(   
                (i+1, PS_point[input_mapping[i]]) for i in range(n_initial+n_final) )
        
        # Retrieve some possibly relevant model parameters
        alpha_s = self.model.get('parameter_dict')['aS']
        mu_r = self.model.get('parameter_dict')['MU_R']   
        
        # Now compute the reduced quantities which will be necessary for evalauting the
        # integrated current
        reduced_PS = counterterm.get_reduced_kinematics(mapped_PS_point)
        
#        misc.sprint('input_mapping %d=%s'%(i_mapping+1, input_mapping))
#        misc.sprint('mapped_flavors=%s'%str(mapped_flavors))
#        misc.sprint('mapped_PS_point=\n%s'%str(mapped_PS_point))
#        misc.sprint('reduced_PS=\n%s'%str(reduced_PS))
#        misc.sprint('counterterm.momenta_dict = %s'%str(counterterm.momenta_dict))
#        misc.sprint('counterterm.process.initial_eg_numbers = %s'%str([l.get('number') for l in counterterm.process.get_initial_legs()]))
#        misc.sprint('counterterm.process.final_numbers = %s'%str([l.get('number') for l in counterterm.process.get_final_legs()]))
        
        # For now the integrated counterterm are always composed of a *single* integrated
        # current, which we can call directly.
        assert(len(counterterm.nodes)==1)
        assert(len(counterterm.nodes[0].nodes)==0)
        integrated_current  = counterterm.nodes[0].current
        reduced_process     = counterterm.process
        momenta_map         = counterterm.momenta_dict
        
        # /!\ Warnings the flavor of the reduced process and well as the current
        # do not match the particular selection. This should however be irrelevant for the
        # evaluation of the counterterm.
        CT_evaluation, all_results = self.all_MEAccessors(
            integrated_current, reduced_PS, reduced_process = reduced_process,
            leg_numbers_map = momenta_map,
            hel_config = hel_config,
            compute_poles = compute_poles )

        all_necessary_ME_calls = []
        for ((spin_index, color_index), current_wgt) in CT_evaluation['values'].items():
            # Now combine the correlators necessary for this current, with those already
            # specified in 'all_necessary_ME_calls'
            all_necessary_ME_calls.append( ( CT_evaluation['spin_correlations'][spin_index],
                                   CT_evaluation['color_correlations'][color_index],
                                   base_objects.EpsilonExpansion(current_wgt) ) )

        # Finally treat the call to the reduced matrix element
        final_weight = base_objects.EpsilonExpansion()
        for (spin_correlators, color_correlators, current_weight) in all_necessary_ME_calls:
            try:
                ME_evaluation, all_ME_results = self.all_MEAccessors(
                   reduced_process, reduced_PS, alpha_s, mu_r,
                   # Let's worry about the squared orders later, we will probably directly fish
                   # them out from the ME_process, since they should be set to a unique combination
                   # at this stage.
                   squared_orders    = None,
                   color_correlation = color_correlators,
                   spin_correlation  = spin_correlators, 
                   hel_config        = hel_config
                )
                #misc.sprint(str(reduced_PS),str(counterterm), reduced_process.nice_string(), spin_correlators, color_correlators, current_weight, ME_evaluation)
            except MadGraph5Error as e:
                logger.critical("""
A reduced matrix element is missing in the library of automatically generated matrix elements.
This is typically what can happen when your process definition is not inclusive over all IR sensitive particles.
Make sure that your process definition is specified using the relevant multiparticle labels (typically 'p' and 'j').
Also make sure that there is no coupling order specification which receives corrections.
The missing process is: %s"""%reduced_process.nice_string())
                raise e
            # Again, for the integrated subtraction counterterms, some care will be needed here
            # for the real-virtual, depending on how we want to combine the two Laurent series.
            final_weight += current_weight*base_objects.EpsilonExpansion(ME_evaluation)

        # Now finally handle the overall prefactor of the counterterm and the counterterm
        # multiplicity factor
        final_weight *= counterterm.prefactor*integrated_CT_multiplicity
            
        # Make sure we have an Laurent series with a double pole at most.
        if not final_weight.max_eps_power()<=0 or \
           not final_weight.min_eps_power()>=-2:
            raise MadGraph5Error("The integrated counterterm:\n%s\n"%str(counterterm)+
                "with reduced matrix element:\n%s\n"%str(reduced_process)+
                "yielded the following laurent series:\n%s\n"%str(final_weight)+
                "which has powers of epsilon not in [-2,-1,0]. This is incorrect for the virtual contribution")   
        
        if compute_poles:
            result = final_weight
        else:
            result = final_weight[0]
        
        # For now always return the original flavor selection, but this may become different
        # when doing initial-final
        return result, flavors
        
    def test_IR_poles(self, test_options):
        """ Compare the IR poles residues in dimensional regularization from the virtual
        contribution and from the integrated counterterm. """
        
        if test_options['seed']:
            random.seed(test_options['seed'])
        
        # Retrieve some possibly relevant model parameters
        alpha_s = self.model.get('parameter_dict')['aS']
        mu_r = self.model.get('parameter_dict')['MU_R']    
        
        # First generate a kinematic point
        # Specifying None forces to use uniformly random generating variables.
        a_virtual_PS_point, _, _, _ = self.phase_space_generator.get_PS_point(None)
        a_virtual_PS_point = phase_space_generators.LorentzVectorDict(
                                   (i+1, mom) for i, mom in enumerate(a_virtual_PS_point) )

        # Now keep track of the results from each process and limit checked
        all_evaluations = {}
        for process_key, (defining_process, mapped_processes) in self.processes_map.items():
            # Make sure that the selected process satisfies the selected process
            if not self.is_part_of_process_selection(
                [defining_process,]+mapped_processes, selection = test_options['process'] ):
                continue
            
            misc.sprint('\nTesting virtual %s'%
                              defining_process.nice_string().replace('Process','process'))

            # Here we use correction_order to select CT subset
            counterterms_to_consider = [
                ct for ct in self.integrated_counterterms[process_key]
                if ct['integrated_counterterm'].count_unresolved() <= 
                                            test_options['correction_order'].count('N') ]
            
            virtual_ME_evaluation, all_results = self.all_MEAccessors(
               defining_process, a_virtual_PS_point, alpha_s, mu_r,
               squared_orders    = None,
               color_correlation = None,
               spin_correlation  = None, 
               hel_config        = None )

            # Turn the evaluation into an epsilon expansion
            virtual_ME_expansion = base_objects.EpsilonExpansion(virtual_ME_evaluation)
            
            # Add the corresponding multiplicity factor if all flavors are to be considered
            if test_options['include_all_flavors']:
                virtual_ME_expansion *= float(1+len(mapped_processes))
            
            # Now loop over all counterterms
            all_integrated_CT_summed_res = base_objects.EpsilonExpansion()
            for i_proc, process in enumerate([defining_process,]+mapped_processes):
                if i_proc > 0 and not test_options['include_all_flavors']:
                    continue
                flavors = process.get_cached_initial_final_pdgs()
#                misc.sprint('Integrated counterterms for flavors: %s'%str(flavors))
                for counterterm in counterterms_to_consider:
#                    misc.sprint('>>Investigating counterterm: %s'%str(counterterm['integrated_counterterm']))
#                    misc.sprint('>>with mappings: %s'%str(counterterm['input_mappings']))
                    for i_mapping, input_mapping in enumerate(counterterm['input_mappings']):
                        # Evaluate the counterterm
                        integrated_CT_res, reduced_flavors = self.evaluate_integrated_counterterm(
                            counterterm, a_virtual_PS_point, flavors, input_mapping, 
                                                          hel_config=None, compute_poles=True)
                        if integrated_CT_res is None:
                            continue
                        misc.sprint('%-20s @ %-35s (permutation #%-2d) => %s' % (
                            str(counterterm['integrated_counterterm']), 
                            str(counterterm['resolved_flavors_combinations'].keys()[0]),
                            i_mapping+1, integrated_CT_res.__str__(format='.16e') ) )
#                        misc.sprint('   :: %s'%(' & '.join('%s => %s'%(str(k),str(v)) for k,v in 
#                                       counterterm['resolved_flavors_combinations'].items() )))
#                        misc.sprint('   :: %s'%(input_mapping))
#                        misc.sprint('   :: %s'%(counterterm['integrated_counterterm'].process.nice_string()))
#                        misc.sprint('   :: leg numbers = %s > %s'%(
#                            ' '.join('%d'%l.get('number') for l in counterterm['integrated_counterterm'].process.get_initial_legs()),
#                            ' '.join('%d'%l.get('number') for l in counterterm['integrated_counterterm'].process.get_final_legs())
#                            ))
                        
                        all_integrated_CT_summed_res += integrated_CT_res
                
            # Add evaluations to the list so as to study how the approximated reals converge towards the real
            evaluation = {
                 'virtual_ME'           : virtual_ME_expansion, 
                 'integrated_CTs'       : all_integrated_CT_summed_res,
                 'defining_process'     : defining_process,
                 'PS_point'             : a_virtual_PS_point }
            
            relative_diff = virtual_ME_expansion.relative_diff(all_integrated_CT_summed_res*-1.)
            # Finite parts are of course expected to differ, so let's not show them
            relative_diff.truncate(min_power = -2, max_power = -1)
            # To be commented out when we will have a full-fledged analysis coded up 
            # in analyze_IR_poles_check()
            misc.sprint('Summary for that PS point:\n%-20s : %s\n%-20s : %s\n%-20s : %s'%(
                 'virtual contrib.',
                 virtual_ME_expansion.__str__(format='.16e'),
                 'integrated CTs.',
                 all_integrated_CT_summed_res.__str__(format='.16e'),
                 'relative diff.',
                 relative_diff.__str__(format='.16e')
            ))
            all_evaluations[process_key] = evaluation
    
        # Now produce a nice output of the evaluations and assess whether this test passed or not.
        return self.analyze_IR_poles_check(all_evaluations, test_options['acceptance_threshold'])    

    def analyze_IR_poles_check(self, all_evaluations, acceptance_threshold):
        """ Analyze the results of the check_IR_pole_residues command. """
        
        #TODO
#        misc.sprint("----- SUMMARY -----")
#        for key, evaluation in all_evaluations.items():
#            misc.sprint("Result for test: %s | %s"%(str(dict(key[0])['PDGs']),key[1]))
#            virtual_ME_expansion = evaluation['virtual_ME']
#            all_integrated_CT_summed_res = evaluation['integrated_CTs']
#            relative_diff = virtual_ME_expansion.relative_diff(all_integrated_CT_summed_res)
#            misc.sprint('\n%-30s : %s\n%-30s : %s\n%-30s : %s'%(
#                 'virtual contrib.',
#                 virtual_ME_expansion.__str__(format='.16e'),
#                 'integrated CTs.',
#                 all_integrated_CT_summed_res.__str__(format='.16e'),
#                 'relative diff.',
#                 relative_diff.__str__(format='.16e')
#            ))
#            This is to be compared with the acceptance_threshold
#            difference_norm = relative_diff.norm()
        
        return True

    def sigma(self, PS_point, process_key, process, flavors, process_wgt, mu_r, mu_f1, mu_f2, *args, **opts):
        """ Overloading of the sigma function from ME7Integrand to include necessary additional contributions. """
        
        sigma_wgt = super(ME7Integrand_V, self).sigma(
                PS_point, process_key, process, flavors, process_wgt, mu_r, mu_f1, mu_f2, *args, **opts)
        
        # This will group all CT results with the same reduced flavors, so
        # as to call the generation-level cuts and observables only once for each
        # flavour configuration.
        # The format is:
        #     { <reduced_flavors_tuple> : [(counterterm, CT_wgt), ] }
        CT_results = {}
        for counterterm_characteristics in self.integrated_counterterms[process_key]:
            for input_mapping in counterterm_characteristics['input_mappings']:
                CT_wgt, reduced_flavors = self.evaluate_integrated_counterterm(
                   counterterm_characteristics, PS_point, flavors, input_mapping,
                   hel_config=None, compute_poles=False)
                    
                if CT_wgt is None or CT_wgt == 0.:
                    continue
    
                # Register this CT in the dictionary CT_results gathering all evaluations
                key = reduced_flavors
                if key not in CT_results:
                    CT_results[key] = [(counterterm_characteristics, CT_wgt)]
                else:
                    CT_results[key].append((counterterm_characteristics, CT_wgt))

        # Now investigate all the different counterterm contributions to add:
        for reduced_flavors, counterterms_characteristics in CT_results.items():
            if not self.pass_flavor_sensitive_cuts(PS_point, reduced_flavors):
                continue
            this_CT_group_wgt = sum(contrib[1] for contrib in counterterms_characteristics)
            if self.apply_observables:
                data_for_observables = {
                    'PS_point'     : PS_point,
                    'flavors'      : reduced_flavors,
                    'counterterms' : counterterms_characteristics }
                self.observable_list.apply_observables(
                                   this_CT_group_wgt*process_wgt, data_for_observables)
                
                # Register this CT_wgt in the global weight.
                sigma_wgt += this_CT_group_wgt

        return sigma_wgt

class ME7Integrand_R(ME7Integrand):
    """ME7Integrand for the computation of a single real-emission type of contribution."""

    divide_by_jacobian = True
    
    def __init__(self, *args, **opts):
        """Initialize a real-emission type of integrand,
        adding additional relevant attributes.
        """
        
        # Initialize the (counter)terms that make up this integrand
        requires = "Constructor of class ME7Integrand_R requires the option '%s'."
        try:
            self.counterterms = opts.pop('counterterms')
        except KeyError:
            raise MadEvent7Error(requires % 'counterterms')
        # Initialize a mapping walker to handle the limits of this integrand
        try:
            self.subtraction_mappings_scheme = opts.pop('subtraction_mappings_scheme')
        except KeyError:
            raise MadEvent7Error(requires % 'subtraction_mappings_scheme')
        try:
            self.walker = mappings.VirtualWalker(self.subtraction_mappings_scheme)
        except KeyError:
            raise MadEvent7Error(
                "Invalid subtraction_mappings_scheme '%s'." %
                self.subtraction_mappings_scheme)
        # Initialize the ME7Integrand
        super(ME7Integrand_R, self).__init__(*args, **opts)
        # Update the initialization inputs
        self.initialization_inputs['options']['counterterms'] = \
            self.counterterms
        self.initialization_inputs['options']['subtraction_mappings_scheme'] = \
            self.subtraction_mappings_scheme

    def get_additional_nice_string_printout_lines(self, format=0):
        """ Return additional information lines for the function nice_string of this integrand."""
        res = []
        if self.counterterms:
            res.append('%-30s:   %d'%('Number of local counterterms', 
               len([1 for CT in sum(self.counterterms.values(),[]) if CT.is_singular()]) ))
        return res
        
    def get_nice_string_process_line(self, process_key, defining_process, format=0):
        """ Return a nicely formated process line for the function nice_string of this 
        integrand."""
        
        GREEN = '\033[92m'
        ENDC = '\033[0m'
        process_string = defining_process.nice_string(print_weighted=False)
        res = GREEN + '  ' + process_string.replace('Process: ','') + ENDC

        if not self.counterterms:
            return res                                                               

        if format < 2:
            if process_key in self.counterterms:
                res += ' | %d local counterterms'%len([
                    1 for CT in self.counterterms[process_key] if CT.is_singular() ])
            else:
                res += ' | 0 local counterterm'
                
        else:
            long_res = [' | with the following local counterterms:']
            for CT in self.counterterms[process_key]:
                if CT.is_singular():
                    if format==2:
                        long_res.append( '   | %s'%CT.__str__(
                            print_n=True, print_pdg=False, print_state=False )  )
                    elif format==3:
                        long_res.append( '   | %s'%CT.__str__(
                            print_n=True, print_pdg=True, print_state=True )  )
                    elif format>3:
                        long_res.append( '   | %s'%str(CT))
            res += '\n'.join(long_res)

        return res

    def combine_color_correlators(self, color_correlators):
        """ This function takes a list of color_correlators in argument, each specified as:
               color_correlator = ( connection_left, connection_right )
            where connection_<x> is given as:
               connection = ( (motherA, emittedA, daughterA), (motherB, emittedB, daughterB), etc...)
            and returns two lists: 
                combined_color_correlators: the list of of new color correlator specifiers 
                                            that arose from the combination of all of them
                multiplier_factors: a list of float factors to multiply the contribution
                                    of each of the combined color correlator returned above.
        """

        # Trivial combination if there is a single one:
        if len(color_correlators)==1:
            return [color_correlators[0],], [1.0,]

        # First weed out the trivial color correlators set to None.
        non_trivial_color_correlators = [ccs for ccs in color_correlators if ccs is not None]
        
        # We don't support terms that factorize a *sum* of color correlator
        if any(len(ccs)!=1 for ccs in non_trivial_color_correlators):
            raise NotImplementedError(
""" The combination of color correlators only the supports the case where each term returned by
    the current factorizes a *single* correlator not a sum of several ones. So re-implement as follows:
    from:
        evaluation = {
              'spin_correlations'  = [ None, ]
              'color_correlations' = [ (correlatorA, correlatorB) ]
              'values'             = [ (0,0) : my_weight ]
        }
    into:
        evaluation = {
              'spin_correlations'  = [ None, ]
              'color_correlations' = [ (correlatorA,), (correlatorB,) ]
              'values'             = [ (0,0) : my_weight,
                                       (0,1) : my_weight ]
        }
""")
        non_trivial_color_correlators = [ccs[0] for ccs in non_trivial_color_correlators]

        # Convert the NLO short-hand convention if present to the general one
        normalized_non_trivial_color_correlators = []
        for cc in non_trivial_color_correlators:
            normalized_non_trivial_color_correlators.append(
                ( cc[0] if not isinstance(cc[0],int) else ((cc[0],-1,cc[0]),),
                  cc[1] if not isinstance(cc[1],int) else ((cc[1],-1,cc[1]),)  )
            )
        non_trivial_color_correlators = normalized_non_trivial_color_correlators
        
        # Place some limitation of the combinations we understand and can support now
        # Only combinations of two color correlators for now
        if len(non_trivial_color_correlators)>2:
            raise NotImplementedError("Only combinations of at most TWO color correlators are supported for now.")
        # Only NLO-type of correlators for now
        if any( len(cc[0])>1 for cc in non_trivial_color_correlators):
            raise NotImplementedError("Only combinations of at most NLO-type of color correlators are supported for now.")
        
        # Treat the easy case first where there is only one or zero color correlator
        if len(non_trivial_color_correlators)==0:
            return [ None,], [1.0,]
        elif len(non_trivial_color_correlators)==1:
            return [ (non_trivial_color_correlators[0],) ], [1.0,]
        
        # Now treat the non-trivial case of the combination of two NLO-type of color correlators
        # For details of this implementation, I refer the reader to the corresponding documentation
        # produced on this topic
        assert(len(non_trivial_color_correlators)==2)
        ccA = non_trivial_color_correlators[0]
        ccB = non_trivial_color_correlators[1]
        # Identify the terms of the compound soft operator S^{i1,j1} \otimes S^{i2,j2} 
        # Each correlator ccA/B will be given in the form of:
        #    (((i, -1, i),), ((j, -1, j),))
        i1_, j1_ = ccA[0][0][0], ccA[1][0][0]
        i2_, j2_ = ccB[0][0][0], ccB[1][0][0]
        
        # Below is my original way of combining the color structure of the two single soft
        def VH_two_iterated_soft_combination((i1,j1),(i2,j2)):
            # now construct all four ways in which the two gluons could have been connected
            # between these four lines.
            if i1==i2:
                connections_A = [
                    ( (i1,-1,i1),(i2,-2,i2) ),
                    ( (i2,-2,i2),(i1,-1,i1) )
                ]
            else:
                connections_A = [
                    ( (i1,-1,i1),(i2,-2,i2) ) if i1>i2 else
                    ( (i2,-2,i2),(i1,-1,i1) )
                ]
            if j1==j2:
                connections_B = [
                    ( (j1,-1,j1),(j2,-2,j2) ),
                    ( (j2,-2,j2),(j1,-1,j1) )
                ]
            else:
                connections_B = [
                    ( (j1,-1,j1),(j2,-2,j2) ) if j1>j2 else
                    ( (j2,-2,j2),(j1,-1,j1) )
                ]
            combined_color_correlators = []
            for connection_A in connections_A:
                for connection_B in connections_B:
                    combined_color_correlators.append(
                        ( connection_A, connection_B )
                    )
            
            multiplier = 1.0/float(len(combined_color_correlators))
            
             # Implementing them as separate calls to the ME, as done in the commented line below,
            # return [ tuple(combined_color_correlators), ], [ multiplier, ]
            # would only yield the same result once the accessor correctly sums over the
            # different color-correlated MEs, so avoid for now.
            
            # Now return the combined color correlators identified
            return [ (ccc,) for ccc in combined_color_correlators ],\
                   [ multiplier, ]*len(combined_color_correlators)
       
        def double_correlator((i,j),(k,l)):
            """ Returns the double correlator of Catani-Grazzini (Eq.113 of hep-ph/9908523v1)
                <M| ( T^-1_i \dot T^-1_j ) * ( T^-1_k \dot T^-1_l ) | M > 
                
            converted into MadGraph's conventions:
            
              ( (a,-1,a),(b,-2,b) ) , ( (c,-1,c),(d,-2,d) ) --> T^-1_a T^-2_b T^-1_c T^-2_d
            """

            # It is important to never commute two color operators acting on the same index, so we must chose carefully which
            # index to pick to carry the gluon index '-2' of the first connection. This can be either 'k' or 'l'.
            if j!=k and j!=l:
                # If all indices are different, we can pick either k or l, it is irrelevant
                index1, index2, index3, index4 = i, k, j, l
            elif j==k and j!=l:
                # If j is equal to k, we must pick l
                index1, index2, index3, index4 = i, l, j, k
            elif j==l and j!=k:
                # If j is equal to l, we must pick k
                index1, index2, index3, index4 = i, k, j, l
            elif j==l and j==k:
                # If j is equal to both l and k, then agin it doesn't matter and we can pick k
                index1, index2, index3, index4 = i, k, j, l
    
            # The sorting according to the first index of each tuple of each of the two convention is to match
            # Madgraph's convention for sorting color connection in the color correlators definition
            return (
                tuple(sorted([ (index1,-1,index1), (index2,-2,index2) ], key = lambda el: el[0], reverse=True)),
                tuple(sorted([ (index3,-1,index3), (index4,-2,index4) ], key = lambda el: el[0], reverse=True))
            )
       
        # Below is a purely abelian combination of the two single softs:
        #  -> second line of Eq. (6.13) of Gabor's hep-ph/0502226v2
        
        def abelian_combination((i1,j1),(i2,j2)):
            # The two terms of the symmetrised sum
            correlator_A   = double_correlator((i1,j1),(i2,j2))
            correlator_B   = double_correlator((i2,j2),(i1,j1))
            # There is an extra factor 2 because in the IR subtraction module, only one combination
            # S(r) S(s) is considered and not the symmetric version S(s) S(r)
            overall_factor = (1.0/4.0)*2.0

            # Group them if equal:
            if correlator_A==correlator_B:
                return [ (correlator_A,), ], [ 2.0*overall_factor ]
            else:
                # Implementing them as separate calls to the ME, as done in the commented line below,
                # return [ (correlator_A,),  (correlator_B,) ], [ overall_factor, overall_factor ]
                # would only yield the same result once the accessor correctly sums over the
                # different color-correlated MEs, so avoid for now.
                return [ (correlator_A,), (correlator_B,) ], [ overall_factor, overall_factor ]

        #return VH_two_iterated_soft_combination((i1_,j1_),(i2_,j2_))
        return abelian_combination((i1_,j1_),(i2_,j2_))

    def combine_spin_correlators(self, spin_correlators):
        """ This function takes several spin-correlators specified int the form
              spin_correlator = ( (legIndex, ( vector_A, vector_B, ...) ),
                                 (anotherLegIndex, ( vector_C, vector_D, ...) ),
                                 etc... ) 
            and returns the list of new correlator specifiers that arises from the
            combination of all of them.
        """
        
        # Trivial combination if there is a single one:
        if len(spin_correlators):
            return spin_correlators[0]

        # This combination is done with a simple concatenation of the lists. Example:
        # Let's say correlator A only correlates to leg #1 with two four-vectors v_A and v_B
        # (these can be seen as a replacement of polarization vectors to consider and summed over)
        #     correlators_A[0]  = ( (1, (v_A, v_B)), )
        # Now correlator B correlates with the two legs #4 and #7 (they must be different 
        #  by construction!), each defined with a single vector v_C and v_B
        #     correlators_B[0]  = ( (4, (v_C,)), (7, (v_D,)) )
        # Then it is clear tha the combined spin correlation should be:
        #     combined_spin_correlator = ( (1, (v_A, v_B)), (4, (v_C,)), (7, (v_D,)) )
        # Notice that this implies that both combinations :
        #    pol_vec_1 = v_A, pol_vec_4 = v_C,  pol_vec_7 = v_D
        # as well as:
        #    pol_vec_1 = v_B, pol_vec_4 = v_C,  pol_vec_7 = v_D
        # will be computed and summed in the resulting spin-correlated matrix element call.
        
        # Make sure the spin correlators don't share common legs
        for i, sc_a in enumerate(spin_correlators):
            for sc_b in spin_correlators[i+1:]:
                assert (len( set(c[0] for c in sc_a)&set(c[0] for c in sc_b))==0 )
        
        return tuple(sum([ (list(sc) if not sc is None else []) for sc in spin_correlators],[]))

    @staticmethod
    def update_all_necessary_ME_calls(all_necessary_ME_calls, new_evaluation):

        new_all_necessary_ME_calls = []
        for ((spin_index, color_index), current_wgt) in new_evaluation['values'].items():
            # Now combine the correlators necessary for this current
            # with those already specified in 'all_necessary_ME_calls'
            for ME_call in all_necessary_ME_calls:
                new_all_necessary_ME_calls.append((
                    # Append this spin correlation to those already present for that call
                    ME_call[0] + [new_evaluation['spin_correlations'][spin_index], ],
                    # Append this color correlation to those already present for that call
                    ME_call[1] + [new_evaluation['color_correlations'][color_index], ],
                    # Append this weight to those already present for that call
                    ME_call[2] + [current_wgt['finite'], ],
                ))
        # Return the new list of necessary ME calls
        return new_all_necessary_ME_calls

    def evaluate_counterterm(
        self, counterterm, PS_point,
        hel_config=None, defining_flavors=None,
        apply_flavour_blind_cuts=True, apply_flavour_cuts=True ):
        """Evaluate a counterterm for a given PS point."""

        # Retrieve some possibly relevant model parameters
        alpha_s = self.model.get('parameter_dict')['aS']
        mu_r = self.model.get('parameter_dict')['MU_R']

        # Now call the walker which hikes through the counterterm structure
        # to return the list of currents and PS points for their evaluation
        # hike_output = {
        #         'currents' : [stroll_output1, stroll_output2, ...],
        #         'matrix_element': (ME_process, ME_PS),
        #         'kinematic_variables' : kinematic_variables (a dictionary) }
        hike_output = self.walker.walk_to_lower_multiplicity(
            PS_point, counterterm, compute_jacobian=self.divide_by_jacobian )
       
        # Access the matrix element characteristics
        ME_process, ME_PS = hike_output['matrix_element']
        
        # Generate what is the kinematics (reduced_PS) returned as a list
        # and the reduced_flavors for this counterterm,
        # by using the defining selected flavors of the real-emission
        # and the real-emission kinematics dictionary
        reduced_PS, reduced_flavors = counterterm.get_reduced_quantities(
            ME_PS, defining_flavors=defining_flavors)

        n_unresolved_left = self.contribution_definition.n_unresolved_particles
        n_unresolved_left -= counterterm.count_unresolved()
        # Apply cuts if requested and return immediately if they do not pass
        if apply_flavour_blind_cuts and not self.pass_flavor_blind_cuts(
            reduced_PS, reduced_flavors,
            n_jets_allowed_to_be_clustered=n_unresolved_left):
            return 0.0, reduced_PS, reduced_flavors
        if apply_flavour_cuts and not self.pass_flavor_sensitive_cuts(
            reduced_PS, reduced_flavors):
            return 0.0, reduced_PS, reduced_flavors

        # The above "hike" can be used to evaluate the currents first and the ME last.
        # Note that the code below can become more complicated when tracking helicities,
        # but let's forget this for now.
        weight = 1.
        assert ((hel_config is None))
        # all_necessary_ME_calls is a list of tuples of the following form:
        # (spin_correlators_to_combine, color_correlators_to_combine, weights_to_combine)
        all_necessary_ME_calls = [ ([], [], []) ]

        for stroll_output in hike_output['currents']:
            stroll_currents = stroll_output['stroll_currents']
            higher_PS_point = stroll_output['higher_PS_point']
            lower_PS_point  = stroll_output['lower_PS_point']
            stroll_vars     = stroll_output['stroll_vars']
            if stroll_vars.has_key('jacobian'):
                weight /= stroll_vars['jacobian']
            for current in stroll_currents:
                # WARNING The use of reduced_process here is fishy (for all but the last)
                current_evaluation, all_current_results = self.all_MEAccessors(
                    current,
                    higher_PS_point=higher_PS_point, lower_PS_point=lower_PS_point,
                    leg_numbers_map=counterterm.momenta_dict,
                    reduced_process=ME_process, hel_config=None,
                    **stroll_vars )
                # Now loop over all spin- and color- correlators required for this current
                # and update the necessary calls to the ME
                if not current['resolve_mother_spin_and_color']:
                    # Make sure no spin- or color-correlations were produced by the current
                    assert(current_evaluation['spin_correlations']==[None,])
                    assert(current_evaluation['color_correlations']==[None,])
                    assert(current_evaluation['values'].keys()==[(0,0),])
                    # WARNING:: this can only work for local 4D subtraction counterterms!
                    # For the integrated ones it is very likely that we cannot use a nested structure,
                    # and there will be only one level anyway,
                    # so there is not need of fancy combination of Laurent series.
                    weight *= current_evaluation['values'][(0,0)]['finite']
                else:
                    all_necessary_ME_calls = ME7Integrand_R.update_all_necessary_ME_calls(
                        all_necessary_ME_calls, current_evaluation)

        # Now perform the combination of the list of spin and color correlators to be merged
        # for each necessary ME call identified
        new_all_necessary_ME_calls = []
        for spin_correlators_to_combine, color_correlators_to_combine, weights_to_combine in all_necessary_ME_calls:
            combined_spin_correlator  = self.combine_spin_correlators(spin_correlators_to_combine)
            # The combination of color correlators can give rise to several ones, each with its multiplier
            combined_color_correlators, multipliers = self.combine_color_correlators(color_correlators_to_combine)
            combined_weight = reduce(lambda x,y: x*y, weights_to_combine)
            for combined_color_correlator, multiplier in zip(combined_color_correlators,multipliers):
                # Finally add the processed combination as a new ME call
                new_all_necessary_ME_calls.append( (
                    combined_spin_correlator,
                    combined_color_correlator,
                    combined_weight*multiplier
                ) )
        all_necessary_ME_calls = new_all_necessary_ME_calls

        # Finally treat the call to the reduced connected matrix elements
        final_weight = 0.0
#        misc.sprint('I got for %s:'%str(counterterm.nice_string()))
        for (spin_correlators, color_correlators, current_weight) in all_necessary_ME_calls:
#            misc.sprint(ME_PS,ME_process.nice_string())
#            misc.sprint(spin_correlators)
#            misc.sprint(color_correlators, current_weight)
            try:
                ME_evaluation, all_ME_results = self.all_MEAccessors(
                   ME_process, ME_PS, alpha_s, mu_r,
                   # Let's worry about the squared orders later, we will probably directly fish
                   # them out from the ME_process, since they should be set to a unique combination
                   # at this stage.
                   squared_orders    = None,
                   color_correlation = tuple(color_correlators) if color_correlators else None,
                   spin_correlation  = tuple(spin_correlators) if spin_correlators else None, 
                   hel_config        = None 
                )
            except MadGraph5Error as e:
                logger.critical("""
A reduced matrix element is missing in the library of automatically generated matrix elements.
This is typically what can happen when your process definition is not inclusive over all IR sensitive particles.
Make sure that your process definition is specified using the relevant multiparticle labels (typically 'p' and 'j').
Also make sure that there is no coupling order specification which receives corrections.
The missing process is: %s"""%ME_process.nice_string())
                raise e
            # for i in ME_PS.keys():
            #     if i in (1, 2): continue
            #     for j in ME_PS.keys():
            #         if j in (1, 2) or j <= i: continue
            #         misc.sprint(i, j, (ME_PS[i]+ME_PS[j]).square())
            # misc.sprint(current_weight, ME_evaluation['finite'])
            # misc.sprint('reduced process = %s' % (
            #     ' '.join('%d(%d)' % (l.get('number'), l.get('id')) for l in
            #              counterterm.process.get_initial_legs()) + ' > ' +
            #     ' '.join('%d(%d)' % (l.get('number'), l.get('id')) for l in
            #              counterterm.process.get_final_legs())
            # ))
            # misc.sprint(counterterm.prefactor)
            # misc.sprint(
            #     'color corr. = %-20s | current = %-20.16f | ME = %-20.16f | Prefactor = %-3f  |  Final = %-20.16f ' % (
            #         str(color_correlators),
            #         current_weight,
            #         ME_evaluation['finite'],
            #         counterterm.prefactor,
            #         current_weight * ME_evaluation['finite'] * counterterm.prefactor
            #     ))
            # Again, for the integrated subtraction counterterms, some care will be needed here
            # for the real-virtual, depending on how we want to combine the two Laurent series.
            final_weight += current_weight*ME_evaluation['finite']
            
        # Now finally handle the overall prefactor of the counterterm
        # and the weight from the disconnected currents
        final_weight *= counterterm.prefactor*weight

        # Returns the corresponding weight and the mapped PS_point.
        # Also returns the mapped_process (for calling the observables), which
        # is typically simply a reference to counterterm.current which is an instance of Process.
        # Notice that the flavors in the ME_process might not be accurate for now.
        return final_weight, reduced_PS, reduced_flavors

    def sigma(self, PS_point, process_key, process, flavors, process_wgt, mu_r, mu_f1, mu_f2, *args, **opts):
        """ Implementation of the short-distance cross-section for the real-emission integrand.
        Counterterms will be computed on top of the actual real-emission integrand."""
        
        # Compute the real-emission matrix element weight in the base ME7Integrand class
        # Notice that the observable will be called already there for the resolved kinematics
        sigma_wgt = super(ME7Integrand_R, self).sigma(PS_point, process_key, process, 
                                  flavors, process_wgt, mu_r, mu_f1, mu_f2, *args, **opts )

        # This will group all CT results with the same reduced kinematics and flavors, so
        # as to call the generation-level cuts and observables only once for each
        # configuration.
        # The format is:
        #     { <tupled_version_of_reduced_PS> :
        #       { 'reduced_PS' : <LorentzVectorList equivalent of the tupled version>,
        #         'flavor_contribs': 
        #           { <reduced_flavors_tuple> : [(counterterm, CT_wgt), ] }
        #       }
        #     }
        CT_results = {}
        for counterterm in self.counterterms[process_key]:
            if not counterterm.is_singular():
                continue

            CT_wgt, reduced_PS, reduced_flavors = self.evaluate_counterterm(
                          counterterm, PS_point, hel_config=None, defining_flavors=flavors,
                          apply_flavour_blind_cuts = True, apply_flavour_cuts = True)
            if CT_wgt == 0.:
                continue

            # Register this CT in the dictionary CT_results gathering all evaluations
            reduced_PS.flags.writeable = False
            key = reduced_PS.data
            if key not in CT_results:
                new_entry = {'flavor_contribs': {}, 'reduced_PS' : reduced_PS}
                CT_results[key] = new_entry
            flavor_results = CT_results[key]['flavor_contribs']
                        
            try:
                flavor_results[reduced_flavors].append( (counterterm, CT_wgt), )
            except KeyError:
                flavor_results[reduced_flavors] = [ (counterterm, CT_wgt), ]

        # Now investigate all the different counterterm contributions to add:
        for reduced_PS_tuple, value in CT_results.items():
            reduced_PS      = value['reduced_PS']
            flavor_contribs = value['flavor_contribs'] 
            for flavor_contrib, counterterm_contribs in flavor_contribs.items():
                this_CT_group_wgt = sum(contrib[1] for contrib in counterterm_contribs)
                # Register this CT_wgt in the global weight.
                sigma_wgt += this_CT_group_wgt
                # Apply observables if requested
                if self.apply_observables:
                    data_for_observables = {
                        'PS_point'     : reduced_PS,
                        'flavors'      : flavor_contrib,
                        'counterterms' : counterterm_contribs }
                    self.observable_list.apply_observables(
                        this_CT_group_wgt*process_wgt, data_for_observables)
        
        return sigma_wgt

    def test_IR_limits(self, test_options):
        """Test how well local counterterms approximate a real-emission matrix element."""

        # Retrieve some possibly relevant model parameters
        alpha_s = self.model.get('parameter_dict')['aS']
        mu_r = self.model.get('parameter_dict')['MU_R']
        # Apply the passed options
        seed = test_options.get('seed', None)
        if seed: random.seed(seed)
        walker_name = test_options.get('walker', None)
        if walker_name is None:
            walker = self.walker
        else:
            walker = mappings.VirtualWalker(walker_name)

        # First generate an underlying Born
        # Specifying None forces to use uniformly random generating variables.
        # Make sure to generate a point within the cuts if necessary:
        max_attempts = 10000
        n_attempts   = 0
        while True:
            n_attempts += 1
            if n_attempts > max_attempts:
                break
            real_emission_PS_point, _, _, _ = self.phase_space_generator.get_PS_point(None)
            if test_options['apply_higher_multiplicity_cuts']:
                if not self.pass_flavor_blind_cuts(
                    real_emission_PS_point,
                    self.processes_map.values()[0][0].get_cached_initial_final_pdgs() ):
                    continue
            break
        if n_attempts > max_attempts:
            raise MadEvent7Error(
                "Could not generate a random kinematic configuration that passes " +
                "the flavour blind cuts in less than %d attempts." % max_attempts )
        n_attempts = 0

        # Loop over processes
        all_evaluations = {}
        for process_key, (defining_process, mapped_processes) in self.processes_map.items():
            misc.sprint("Considering %s"%defining_process.nice_string())
            # Make sure that the selected process satisfies the selection requirements
            if not self.is_part_of_process_selection(
                [defining_process, ]+mapped_processes,
                selection=test_options['process'] ):
                continue
            
            a_real_emission_PS_point = copy.copy(real_emission_PS_point)
            while (test_options['apply_higher_multiplicity_cuts'] and
                   not self.pass_flavor_sensitive_cuts(
                       a_real_emission_PS_point,
                       defining_process.get_cached_initial_final_pdgs() ) ):
                n_attempts += 1
                if n_attempts > max_attempts:
                    break
                a_real_emission_PS_point, _, _, _ = self.phase_space_generator.get_PS_point(None)
            if n_attempts > max_attempts:
                raise MadEvent7Error(
                    "Could not generate a random kinematic configuration that passes " +
                    "the flavour blind cuts in less than %d attempts." % max_attempts )
            n_attempts = 0

            # Make sure to have the PS point provided as LorentzVectorDict
            a_real_emission_PS_point = phase_space_generators.LorentzVectorDict(
                (i+1, mom) for i, mom in enumerate(a_real_emission_PS_point) )

            # Use correction_order to select CT subset
            counterterms_to_consider = [
                ct for ct in self.counterterms[process_key]
                if ct.count_unresolved() <= test_options['correction_order'].count('N') ]
            
            # Select the limits to be probed interpreting limits as a regex pattern.
            # If no match is found, then reconstruct the singular structure from the limits
            # provided
            selected_counterterms = self.find_counterterms_matching_regexp(
                                counterterms_to_consider, test_options['limits'] )
            if selected_counterterms:
                selected_singular_structures = [
                    ct.reconstruct_complete_singular_structure()
                    for ct in selected_counterterms]
            else:
                limit_str = test_options['limits']
                ss = mappings.sub.SingularStructure.from_string(
                    limit_str, defining_process)
                if ss is None:
                    logger.critical(
                        "%s is not a valid limits specification" % limit_str)
                    return
                selected_singular_structures = [ss, ]

            misc.sprint('Reconstructed complete singular structure: \n'+'\n'.join(
                str(ss) for ss in selected_singular_structures ))

            # Loop over approached limits
            process_evaluations = {}
            for limit in selected_singular_structures:
                misc.sprint("Approaching limit %s" % str(limit) )
                # Select counterterms to evaluate
                counterterms_to_evaluate = [ct for ct in counterterms_to_consider]
                if test_options['counterterms']:
                    counterterm_pattern = test_options['counterterms']
                    if counterterm_pattern.startswith('def'):
                        counterterm_pattern = str(limit)
                    counterterms_to_evaluate = self.find_counterterms_matching_regexp(
                                                  counterterms_to_evaluate, counterterm_pattern )
                # Progressively approach the limit, using a log scale
                limit_evaluations = {}
                n_steps = test_options['n_steps']
                min_value = test_options['min_scaling_variable']
                base = min_value ** (1./n_steps)
                for step in range(n_steps+1):
                    # Determine the new phase-space point
                    scaling_parameter = base ** step
                    scaled_real_PS_point = walker.approach_limit(
                        a_real_emission_PS_point,
                        limit, scaling_parameter, defining_process )
                    # Initialize result
                    this_eval = {}
                    # Evaluate ME
                    ME_evaluation, all_results = self.all_MEAccessors(
                       defining_process, scaled_real_PS_point, alpha_s, mu_r,
                       squared_orders    = None,
                       color_correlation = None,
                       spin_correlation  = None,
                       hel_config        = None )
                    this_eval['ME'] = ME_evaluation['finite']
                    misc.sprint('For scaling variable %.3e, weight from ME = %.16f' %(
                                              scaling_parameter, ME_evaluation['finite'] ))
                    # Loop over counterterms
                    for counterterm in counterterms_to_evaluate:
                        ct_weight, _, _ = self.evaluate_counterterm(
                            counterterm,
                            scaled_real_PS_point, 
                            hel_config=None,
                            apply_flavour_blind_cuts=test_options['apply_lower_multiplicity_cuts'],
                            apply_flavour_cuts=test_options['apply_lower_multiplicity_cuts'] )
                        this_eval[str(counterterm)] = ct_weight
                        misc.sprint('Weight from CT %s = %.16f' %
                                    (str(counterterm), ct_weight) )
                        misc.sprint('Ratio: %.16f'%( ct_weight/float(ME_evaluation['finite']) ))
                    limit_evaluations[scaling_parameter] = this_eval

                process_evaluations[str(limit)] = limit_evaluations

            process_string = defining_process.base_string()
            if defining_process.has_key('n_loops'):
                process_string += " @ " + str(defining_process['n_loops']) + " loops"
            all_evaluations[process_string] = process_evaluations

        # Now produce a nice matplotlib of the evaluations
        # and assess whether this test passed or not
        return self.analyze_IR_limits_test(
            all_evaluations, test_options['acceptance_threshold'],
            seed=seed, show=test_options['show_plots'], save=test_options['save_plots'] )

    @staticmethod
    def analyze_IR_limit(
        evaluations, acceptance_threshold,
        title=None, def_ct=None, plot_all=True, show=True, filename=None ):

        import matplotlib.pyplot as plt

        test_failed = False

        # Produce a plot of all counterterms
        x_values = sorted(evaluations.keys())
        lines = evaluations[x_values[0]].keys()
        # Skip ME-def line if there is no defining ct
        plot_def = def_ct and def_ct in lines

        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')

        plt.figure(1)
        if title: plt.title(title)
        plt.xlabel('$\lambda$')
        plt.ylabel('Integrands')
        plt.xscale('log')
        plt.yscale('log')
        plt.grid(True)
        total           = [0., ] * len(x_values)
        ME_minus_def_ct = [0., ] * len(x_values)
        for line in lines:
            y_values = [abs(evaluations[x][line]) for x in x_values]
            for i in range(len(x_values)):
                total[i] += evaluations[x_values[i]][line]
            if plot_def and (line == "ME" or line == def_ct):
                def_ct_sign = 1
                if line == def_ct:
                    def_ct_sign = (-1) ** def_ct.count("(")
                for i in range(len(x_values)):
                    ME_minus_def_ct[i] += def_ct_sign * evaluations[x_values[i]][line]
            if plot_all:
                if '(' in line:
                    style = '--'
                else:
                    style = '-'
                plt.plot(x_values, y_values, style, label=line)
        if plot_def:
            abs_ME_minus_def_ct = [abs(y) for y in ME_minus_def_ct]
            plt.plot(x_values, abs_ME_minus_def_ct, label='ME-def')
        abs_total = [abs(y) for y in total]
        plt.plot(x_values, abs_total, label='TOTAL')
        plt.legend()
        if filename:
            plt.savefig(filename + '_integrands.pdf')

        plt.figure(2)
        if title: plt.title(title)
        plt.xlabel('$\lambda$')
        plt.ylabel('Ratio to ME')
        plt.xscale('log')
        plt.grid(True)
        for line in lines:
            y_values = [abs(evaluations[x][line]/evaluations[x]["ME"]) for x in x_values]
            if '(' in line:
                style = '--'
            else:
                style = '-'
            plt.plot(x_values, y_values, style, label=line)
        plt.legend()
        if filename:
            plt.savefig(filename + '_ratios.pdf')

        # Check that the ratio of def_ct to the ME is close to -1
        if plot_def and not test_failed:
            def_ct_2_ME_ratio = evaluations[x_values[0]][def_ct]
            def_ct_2_ME_ratio /= evaluations[x_values[0]]["ME"]
            foo_str = "The ratio of the defining CT to the ME at lambda = %s is: %s."
            print foo_str % (x_values[0], def_ct_2_ME_ratio)
            test_failed = abs(def_ct_2_ME_ratio+1) > acceptance_threshold
        # Check that the ratio between total and ME is close to 0
        if not test_failed:
            total_2_ME_ratio = total[0]
            total_2_ME_ratio /= evaluations[x_values[0]]["ME"]
            foo_str = "The ratio of the total to the ME at lambda = %s is: %s."
            print foo_str % (x_values[0], total_2_ME_ratio)
            test_failed = abs(total_2_ME_ratio) > acceptance_threshold

        plt.figure(3)
        if title: plt.title(title)
        plt.xlabel('$\lambda$')
        plt.ylabel('Weighted integrands')
        plt.xscale('log')
        plt.yscale('log')
        plt.grid(True)
        if plot_def:
            wgt_ME_minus_def_ct = [abs(x_values[i] * ME_minus_def_ct[i])
                                   for i in range(len(x_values))]
            plt.plot(x_values, wgt_ME_minus_def_ct, label='ME-def')
        wgt_total = [abs(x_values[i] * total[i]) for i in range(len(x_values))]
        plt.plot(x_values, wgt_total, label='TOTAL')
        plt.legend()
        if filename:
            plt.savefig(filename + '_weighted.pdf')

        if show:
            plt.show()
        else:
            plt.close('all')

        return not test_failed

    def analyze_IR_limits_test(
        self, all_evaluations, acceptance_threshold,
        seed=None, show=True, save=False):
        """Analyze the results of the test_IR_limits command."""

        test_failed = False
        results = dict()
        for (process, process_evaluations) in all_evaluations.items():
            results[process] = dict()
            for (limit, limit_evaluations) in process_evaluations.items():
                proc, loops = process.split("@")
                title = "$" + proc + "$"
                title = title.replace('~','x').replace('>','\\to').replace(' ','\\;')
                title = title.replace('+','^+').replace('-','^-')
                title += "@" + loops + " approaching " + limit
                if seed: title += " (seed %d)" % seed
                filename = None
                if save:
                    filename = copy.copy(process)
                    filename = filename.replace(">", "_")
                    filename = filename.replace(" ", "").replace("~", "x")
                    filename = filename.replace("@", "_") + "_" + limit
                    if seed: filename += str(seed)
                results[process][limit] = self.analyze_IR_limit(
                    limit_evaluations, acceptance_threshold=acceptance_threshold,
                    title=title, def_ct=limit, show=show, filename=filename )
                if not results[process][limit]:
                    test_failed = True
        if test_failed:
            print "analyse_IR_limits_test result:"
            for process in results.keys():
                print "    " + str(process)
                for limit in results[process].keys():
                    if results[process][limit]:
                        print "    " * 2 + str(limit) + ": passed"
                    else:
                        print "    " * 2 + str(limit) + ": NOT passed"
            return False
        return True
    
class ME7Integrand_RR(ME7Integrand_R):
    """ ME7Integrand for the computation of a double real-emission type of contribution."""
    def sigma(self, PS_point, process_key, process, flavors, process_wgt, mu_r, mu_f1, mu_f2, *args, **opts):
        return super(ME7Integrand_RR, self).sigma(PS_point, process_key, process, flavors, process_wgt, 
                                                                        mu_r, mu_f1, mu_f2, *args, **opts)

class ME7Integrand_RRR(ME7Integrand_R):
    """ ME7Integrand for the computation of a double real-emission type of contribution."""
    def sigma(self, PS_point, process_key, process, flavors, process_wgt, mu_r, mu_f1, mu_f2, *args, **opts):
        return super(ME7Integrand_RRR, self).sigma(PS_point, process_key, process, flavors, process_wgt, 
                                                                        mu_r, mu_f1, mu_f2, *args, **opts)

#===============================================================================
# ME7IntegrandList
#===============================================================================
class ME7IntegrandList(base_objects.PhysicsObjectList):
    """ Container for ME7Integrnds."""
    
        
    integrands_natural_order = [
        ('LO',    (ME7Integrand_B, ME7Integrand_LIB) ),
        ('NLO',   (ME7Integrand_R, ME7Integrand_V) ),
        ('NNLO',  (ME7Integrand_RR, ) ),
        ('NNNLO', (ME7Integrand_RRR, ) )
    ]
    
    def is_valid_element(self, obj):
        """Test if object obj is a valid instance of ME7Integrand."""
        return isinstance(obj, ME7Integrand)
    
    def get_integrands_of_order(self, correction_order):
        """ Returns a list of all contributions of a certain correction_order in argument."""
        return ME7IntegrandList([integrand for integrand in self if
                integrand.contribution_definition.correction_order==correction_order])

    def get_integrands_of_type(self, correction_classes):
        """ Returns a list of all contributions that are direct instances of certain classes."""
        if not isinstance(correction_classes, tuple):
            if isinstance(correction_classes, list):
                correction_classes = tuple(correction_classes)
            else:
                correction_classes = (correction_classes,)                
        return ME7IntegrandList([integrand for integrand in self if 
                                                isinstance(integrand, correction_classes)])

    def nice_string(self, format=0):
        """ A nice representation of a list of contributions. 
        We can reuse the function from ContributionDefinitions."""
        return base_objects.ContributionDefinitionList.contrib_list_string(self, 
                                                                            format=format)

    def sort_integrands(self):
        """ Sort integrands according to the order dictated by the class attribute
        'integrands_natural_order'"""
        
        new_order = []
        for correction_order, integrand_types in self.integrands_natural_order:
            for integrand_type in integrand_types:
                selected_integrands = self.get_integrands_of_order(correction_order).\
                                        get_integrands_of_type(integrand_type)
                new_order.extend(selected_integrands)
                for integrand in selected_integrands:
                    self.pop(self.index(integrand))
    
        # Finally all remaining contributions of unknown types
        new_order.extend(self)
        self[:] = new_order


# Integrand classes map is defined here as module variables. This map can be overwritten
# by the interface when using a PLUGIN system where the user can define his own Integrand.
# Notice that this must be placed after all the Integrand daughter classes in this module have been declared.
ME7Integrand_classes_map = {'Born': ME7Integrand_B,
                            'LoopInduced_Born': ME7Integrand_LIB,
                            'Virtual': ME7Integrand_V,
                            'SingleReals': ME7Integrand_R,
                            'DoubleReals': ME7Integrand_RR,
                            'TripleReals': ME7Integrand_RRR,
                            'Unknown': None}
