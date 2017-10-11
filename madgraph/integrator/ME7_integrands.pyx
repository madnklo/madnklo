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
"""Collection of ME7 Integrand classes to be compiled with Cython
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
# ME7CythonIntegrand
#===============================================================================
class ME7CythonIntegrand(integrands.VirtualIntegrand):
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
        if cls is ME7CythonIntegrand:
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
            target_class = ME7CythonIntegrand_classes_map[target_type]

            if not target_class:
                raise MadGraph5Error("Could not determine the class of integrand of type '%s' to be added for"%target_type+
                                     " the contribution definiton:\n%s"%str(contribution_definition.nice_string()))

            return super(ME7CythonIntegrand, cls).__new__(target_class, *all_args, **opt)
        else:
            return super(ME7CythonIntegrand, cls).__new__(cls, *all_args, **opt)
    
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
        
        super(ME7CythonIntegrand, self).__init__()
        
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
        res = ['%-30s:   %s'%('ME7CythonIntegrand_type',type(self))]
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
            raise MadGraph5Error("The ME7CythonIntegrand must be initialized with a ModelReader instance.")

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

    def find_counterterms_matching_limit_type_with_regexp(self, counterterms, limit_type=None):
        """ Find all mappings that match a particular limit_type given in argument
        (takes a random one if left to None). This function is placed here given that
        it can be useful for both the ME7Integrnd_V and ME7_integrand_R."""

        # First select only the counterterms which are not pure matrix elements 
        # (i.e. they have singular structures).
        selected_counterterms = [ct for ct in counterterms if ct.is_singular()]
        
        if len(selected_counterterms)==0:
            return []

        returned_counterterms = []
        if not limit_type:
            returned_counterterms.append(random.choice(selected_counterterms))
        else:
            for counterterm in selected_counterterms:
                if re.match(limit_type, str(counterterm)):
                    returned_counterterms.append(counterterm)

        return returned_counterterms

    def pass_flavor_blind_cuts(self, input_PS_point, process_pdgs):
        """ Implementation of a minimal set of isolation cuts. This can be made much nicer in the future and 
        will probably be taken outside of this class so as to ease user-specific cuts, fastjet, etc...
        We consider here a two-level cuts system, this first one of which is flavour blind.
        This is of course not IR safe at this stage!"""
        
        # This is a temporary function anyway which should eventually be replaced by a full
        # fledged module for handling generation level cuts, which would make use of fj-core.
        # The PS point in input is sometimes provided as a dictionary or a flat list, but
        # we need it as a flat list here, so we force the conversion
        PS_point = input_PS_point.to_list()

        # These cuts are not allowed to resolve flavour, but only whether a particle is a jet or not
        def is_a_jet(pdg):
            return pdg in range(1,7)+range(-1,-7,-1)+[21]
        
        if __debug__: logger.debug( "Processing flavor-blind cuts for process %s and PS point:\n%s"%(
                        str(process_pdgs), PS_point.__str__(n_initial=self.phase_space_generator.n_initial) ))

        pt_cut = self.run_card['ptj']
        dr_cut = self.run_card['drjj']
                
        if pt_cut <= 0. and dr_cut <= 0.:
            return True

        # Apply the Ptj cut first
        for i, p in enumerate(PS_point[self.n_initial:]):
            if not is_a_jet(process_pdgs[1][i]):
                continue
            if __debug__: logger.debug('p_%i.pt()=%.5e'%((self.n_initial+i),p.pt()))
            if p.pt() < pt_cut:
                return False

        # And then the drjj cut
        for i, p1 in enumerate(PS_point[self.n_initial:]):
            for j, p2 in enumerate(PS_point[self.n_initial+i+1:]):
                if not is_a_jet(process_pdgs[1][i]) and\
                   not is_a_jet(process_pdgs[1][i+1+j]):
                    continue
                if __debug__: logger.debug('deltaR(p_%i,p_%i)=%.5e'%(
                     self.n_initial+i, self.n_initial+i+1+j, p1.deltaR(p2)))
                if p1.deltaR(p2) < dr_cut:
                    return False

        return True

    def pass_flavor_sensitive_cuts(self, input_PS_point, flavors):
        """ Implementation of a minimal set of isolation cuts. This can be made much nicer in the future and 
        will probably be taken outside of this class so as to ease user-specific cuts, fastjet, etc...
        We consider here a two-level cuts system, this second one of which is flavour sensitive."""

        PS_point = input_PS_point.to_list()

        if __debug__: logger.debug( "Processing flavor-sensitive cuts for flavors %s and PS point:\n%s"%(
                        str(flavors), PS_point.__str__(n_initial=self.phase_space_generator.n_initial) ))

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
                
        if __debug__: logger.debug("Considering the following PS point:\n%s"%(PS_point.__str__(
                                            n_initial=self.phase_space_generator.n_initial) ))
        
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
        for process_key, (process, mapped_processes) in self.processes_map.items():
            if __debug__: logger.debug('Now considering the process group from %s.'%process.nice_string())
            
            this_process_wgt = wgt
            process_pdgs = ( tuple(process.get_initial_ids()),
                             tuple(process.get_final_ids_after_decay()) )

            # Apply flavor blind cuts
            if not self.pass_flavor_blind_cuts(PS_point, process_pdgs):
                if __debug__: logger.debug('Event failed the flavour_blind generation-level cuts.')
                if __debug__: logger.debug(misc.bcolors.GREEN + 'Returning a weight of 0. for this integrand evaluation.' + misc.bcolors.ENDC)            
                return 0.0
            
            all_processes = [process,]+mapped_processes
            all_process_pdgs = [ ( tuple(process.get_initial_ids()),
                   tuple(process.get_final_ids_after_decay()) ) for proc in all_processes ]
            # Add mirror processes if present
            all_process_pdgs.extend([ ( tuple(proc.get_initial_ids()[::-1]),
                tuple(proc.get_final_ids()) )  for proc in all_processes if proc.get('has_mirror_process')])

            if abs(self.run_card['lpp1'])==abs(self.run_card['lpp2'])==1:
                if __debug__: logger.debug("Bjorken x's x1, x2, sqrt(x1*x2*s): %.5e, %.5e. %.5e"%(
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

        ME_evaluation, all_results = self.all_MEAccessors(process, PS_point, alpha_s, mu_r, pdgs=flavors)

        sigma_wgt *= ME_evaluation['finite']
        
        if self.apply_observables:
            data_for_observables = {'PS_point': PS_point, 'flavors' : flavors}
            self.observable_list.apply_observables(sigma_wgt*process_wgt, data_for_observables)

        return sigma_wgt

##########################################################################################
# Daughter classes for the various type of ME7Integrands which require a redefinition 
# of sigma() and possibly other functions.
##########################################################################################

##########################################################################################
# ME7CythonIntegrand_B
##########################################################################################
class ME7CythonIntegrand_B(ME7CythonIntegrand):
    """ME7CythonIntegrand for the computation of a Born type of contribution."""
    
    def sigma(self, PS_point, process_key, process, flavors, process_wgt,
                                                        mu_r, mu_f1, mu_f2, *args, **opts):

        return super(ME7CythonIntegrand_B, self).sigma(
            PS_point, process_key, process, flavors, process_wgt,
            mu_r, mu_f1, mu_f2, *args, **opts
        )

##########################################################################################
# ME7CythonIntegrand_LIB
##########################################################################################
class ME7CythonIntegrand_LIB(ME7CythonIntegrand):
    """ME7CythonIntegrand for the computation of a Born type of contribution."""
    
    def sigma(self, PS_point, process_key, process, flavors, process_wgt,
                                                        mu_r, mu_f1, mu_f2, *args, **opts):

        return super(ME7CythonIntegrand_LIB, self).sigma(
            PS_point, process_key, process, flavors, process_wgt,
            mu_r, mu_f1, mu_f2, *args, **opts
        )
        
##########################################################################################
# ME7CythonIntegrand_R
##########################################################################################
class ME7CythonIntegrand_R(ME7CythonIntegrand):
    """ME7CythonIntegrand for the computation of a Born type of contribution."""
    
    def sigma(self, PS_point, process_key, process, flavors, process_wgt,
                                                        mu_r, mu_f1, mu_f2, *args, **opts):

        return super(ME7CythonIntegrand_R, self).sigma(
            PS_point, process_key, process, flavors, process_wgt,
            mu_r, mu_f1, mu_f2, *args, **opts
        )

##########################################################################################
# ME7CythonIntegrand_V
##########################################################################################
class ME7CythonIntegrand_V(ME7CythonIntegrand):
    """ME7CythonIntegrand for the computation of a Born type of contribution."""
    
    def sigma(self, PS_point, process_key, process, flavors, process_wgt,
                                                        mu_r, mu_f1, mu_f2, *args, **opts):

        return super(ME7CythonIntegrand_V, self).sigma(
            PS_point, process_key, process, flavors, process_wgt,
            mu_r, mu_f1, mu_f2, *args, **opts
        )
        
##########################################################################################
# ME7CythonIntegrand_RR
##########################################################################################
class ME7CythonIntegrand_RR(ME7CythonIntegrand):
    """ME7CythonIntegrand for the computation of a Born type of contribution."""
    
    def sigma(self, PS_point, process_key, process, flavors, process_wgt,
                                                        mu_r, mu_f1, mu_f2, *args, **opts):

        return super(ME7CythonIntegrand_RR, self).sigma(
            PS_point, process_key, process, flavors, process_wgt,
            mu_r, mu_f1, mu_f2, *args, **opts
        )
        
##########################################################################################
# ME7CythonIntegrand_RRR
##########################################################################################
class ME7CythonIntegrand_RRR(ME7CythonIntegrand):
    """ME7CythonIntegrand for the computation of a Born type of contribution."""
    
    def sigma(self, PS_point, process_key, process, flavors, process_wgt,
                                                        mu_r, mu_f1, mu_f2, *args, **opts):

        return super(ME7CythonIntegrand_RRR, self).sigma(
            PS_point, process_key, process, flavors, process_wgt,
            mu_r, mu_f1, mu_f2, *args, **opts
        )
        
# Integrand classes map is defined here as module variables. This map can be overwritten
# by the interface when using a PLUGIN system where the user can define his own Integrand.
# Notice that this must be placed after all the Integrand daughter classes in this module have been declared.
ME7CythonIntegrand_classes_map = {'Born': ME7CythonIntegrand_B,
                            'LoopInduced_Born': ME7CythonIntegrand_LIB,
                            'Virtual': ME7CythonIntegrand_V,
                            'SingleReals': ME7CythonIntegrand_R,
                            'DoubleReals': ME7CythonIntegrand_RR,
                            'TripleReals': ME7CythonIntegrand_RRR,
                            'Unknown': None}

