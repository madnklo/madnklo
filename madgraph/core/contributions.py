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
"""Classes for implementions Contributions.
Contributions are abstract layer, between the madgraph_interface and amplitudes.
They typically correspond to the Born, R, V, RV, RR, VV, etc.. pieces of higher
order correction computations.
"""
import array
import copy
import itertools
import logging
import os
pjoin = os.path.join

import madgraph.core.base_objects as base_objects
import madgraph.core.diagram_generation as diagram_generation
import madgraph.iolibs.export_v4 as export_v4
import madgraph.various.misc as misc
from madgraph import InvalidCmd, MadGraph5Error

logger = logging.getLogger('madgraph.contributions')

#===============================================================================
# Contribution mother class
#===============================================================================
class Contribution(object):
    """ An high-level object that stores all the relevant information for a 
    particular contribution to a cross-section calculation, including higher
    order corrections.
    We do not inherit from PhysicsObject as we want to depart from the dict
    structure imposed by this mother."""
    
    def __new__(cls, contribution_definition, MG_options, **opt):
        if cls is Contribution:
            target_class = None
            if contribution_definition.correction_order == 'LO':
                if contribution_definition.n_loops == 0 and \
                  contribution_definition.n_unresolved_particles == 0:
                   target_class = Contribution_B
                elif contribution_definition.n_loops == 0 and \
                  contribution_definition.n_unresolved_particles == 0:
                   target_class = Contribution_LI                    
            elif contribution_definition.correction_order == 'NLO':
                if contribution_definition.n_loops == 1 and \
                   contribution_definition.n_unresolved_particles == 0:
                    target_class = Contribution_V
                elif contribution_definition.n_loops == 0 and \
                     contribution_definition.n_unresolved_particles == 1:
                    target_class = Contribution_R
            elif contribution_definition.correction_order == 'NNLO':
                raise MadGraph5Error("NNLO type of contributions not implemented yet.")                
            if not target_class:
                raise MadGraph5Error("Could not determine the type of contribution to be added for"+
                                     " the contribution definiton:\n%s"%str(contribution_definition))
            return super(Contribution, cls).__new__(target_class, contribution_definition, MG_options, **opt)
        else:
            return super(RunCard, cls).__new__(cls, finput, **opt)
    
    def __init__(self, contribution_definition, MG_options, **opts):
        """ Instantiates a particular contribution."""
        self.contribution_definition = contribution_definition
        self.amplitudes              = diagram_generation.AmplitudeList()
        self.options                 = MG_options
        self._curr_exporter          = None
        # The following two attributes dicate the type of Exporter which will be assigned to this contribution
        self.output_type             = 'default'
        self.format                  = 'standalone'
        self._export_dir             = 'None'
        if self.options['group_subprocesses'] == 'Auto':
            self.collect_mirror_procs = True
        else:
            self.collect_mirror_procs   = self.options['group_subprocesses']
        # Options relevant only for LO diagram generation
        self.ignore_six_quark_processes = False
        self.diagram_filter             = False 
        self.optimize                   = False

    def set_export_dir(self, prefix):
        """ Assigns an export directory name."""
        dir_name = self.contribution_definition.get_shell_name()
        # Use the name of the first process since we can't access the name of the ProcessDefinition
        dir_name += '_%s'%self.amplitudes[0].get('process').shell_string(
            schannel=False, forbid=False, main=False, pdg_order=False, print_id = False)
        dir_name += '_%d'%self.contribution_definition.process_definition.get('id')
        export_dir = pjoin(prefix, dir_name)
        if os.path.exists(export_dir):
            raise MadGraph5Error("The following contribution:\n"+self.nice_string()+
                "\ncannot be exported at location:\n"+self._curr_exporter+
                "\nsince this directory already exists.")
        self._export_dir = export_dir

    def initialize_exporter(self, cmd_interface, noclean, group_subprocesses=True):
        """ Initialize the exporter that will be associated to that particular contribution.
        noclean specifies what to do in case the output directory already exists and group_subprocesses
        whether the exporter should attempt to group identical subprocesses.
        """
        self.set_export_dir(cmd_interface._export_dir)
        self._curr_exporter = export_v4.ExportV4Factory(
                cmd_interface, noclean, output_type=self.output_type, group_subprocesses=group_subprocesses,
                curr_amps = self.amplitudes,
                export_dir = self._export_dir,
                format = self.format)
    
    def copy_template(self, model):
        """ Copy the template structure for that contribution. Quite often, this limits itself to aksing its
        exporter to do this."""
        return self._curr_exporter.copy_template(model)

    def pass_information_from_cmd(self, cmd_interface):
        """ Pass information from the command_interface to this contribution. Most of the time, this only amounts
        to passing information to the active exporter."""
        return self._curr_exporter.pass_information_from_cmd(cmd_interface)
    
    def export(self, nojpeg=False, group_processes=True, args=[]):
        """ Perform the export duties, that include generation of the HelasMatrixElements and 
        the actual output of the matrix element code by the exporter."""
        return
        raise NotImplementedError

    def finalize(self, flaglist=[], interface_history=[]):
        """ Finalize the output of the code necessary for this contribution."""
        return
        raise NotImplementedError
    
    def nice_string(self):
        """ Nice string representation of self."""
        BLUE = '\033[94m'
        GREEN = '\033[92m'
        ENDC = '\033[0m'
        res = ['%-30s:   %s'%('contribution_type',type(self))]
        res.extend([self.contribution_definition.nice_string()])
        if self.amplitudes:
            res.append('Amplitudes generated for the following processes:')
            for amp in self.amplitudes:
                res.append(BLUE+'  %s'%amp.get('process').nice_string().replace('Process: ','')+ENDC)
        else:
            res.append(BLUE+'No amplitudes generated yet.'+ENDC)                
        return '\n'.join(res)
        
    def generate_amplitudes(self):
        """ Generates the relevant amplitudes for this contribution."""

        myproc = diagram_generation.MultiProcess(self.contribution_definition.process_definition,
                    collect_mirror_procs = self.collect_mirror_procs,
                    ignore_six_quark_processes = self.ignore_six_quark_processes,
                    optimize=self.optimize, diagram_filter=self.diagram_filter)

        for amp in myproc.get('amplitudes'):
            if amp not in self.amplitudes:
                self.amplitudes.append(amp)
            else:
                logger.warning('Duplicate process found in contribution '+
                               '%s. Sanity check needed.'%str(type(self)))
 
class Contribution_B(Contribution):
    """ Implements the handling of a Born-type of contribution."""
    
    def __init__(self, contribution_definition, MG_options, **opts):
        super(Contribution_B,self).__init__(contribution_definition, MG_options, **opts)
        self.ignore_six_quark_processes = self.options['ignore_six_quark_processes'] if \
            "ignore_six_quark_processes" in self.options else []
        if 'diagram_filter' in opts:
            self.diagram_filter = opts['diagram_filter']
        if 'optimize' in opts:
            self.optimize = opts['optimize']

class Contribution_LIB(Contribution_B):
    """ Implements the handling of loop-induced Born-type of contribution."""
    pass

class Contribution_R(Contribution):
    """ Implements the handling of a real-emission type of contribution."""
    pass

class Contribution_V(Contribution):
    """ Implements the handling of a virtual type of contribution."""
    pass

class ContributionList(base_objects.PhysicsObjectList):
    """ A container for storing a list of contributions."""
    
    def is_valid_element(self, obj):
        """Test if object obj is a valid instance of Contribution."""
        return isinstance(obj, Contribution)
    
    def get_contributions_of_order(self, correction_order):
        """ Returns a list of all contributions of a certain correction_order in argument."""
        return ContributionList([contrib for contrib in self if
                contrib.contribution_definition.correction_order==correction_order])

    def get_contributions_of_type(self, correction_classes):
        """ Returns a list of all contributions that are direct instances of certain classes."""
        if not isinstance(correction_classes, tuple):
            if isinstance(correction_classes, list):
                correction_classes = tuple(correction_classes)
            else:
                correction_classes = (correction_classes,)                
        return ContributionList([contrib for contrib in self if isinstance(contrib, correction_classes)])
    
    def nice_string(self):
        """ A nice representation of a list of contributions. 
        We can reuse the function from ContributionDefinitions."""
        return base_objects.ContributionDefinitionList.contrib_list_string(self)

    def generate_amplitudes(self):
        """ Generate all amplitudes for the contributions in this list."""
        
        ordered_contributions = [
            ('LO',   (Contribution_B, Contribution_LIB) ),
            ('NLO',  (Contribution_R, Contribution_V) ),
            ('NNLO', () )
        ]
        
        for correction_order, contribution_types in ordered_contributions:
            for contrib_type in contribution_types:
                selected_contribs = self.get_contributions_of_order(correction_order).\
                                        get_contributions_of_type(contrib_type)
                if selected_contribs:
                    logger.info('Generating diagrams for the %d contribution%s of type %s (%s)...'%(
                        len(selected_contribs), 's' if len(selected_contribs)>1 else '',
                        contrib_type.__name__, correction_order))
                for i, contrib in enumerate(selected_contribs):
                    logger.info('%s (%s) %d/%d'%
                        (contrib_type.__name__, correction_order, i+1, len(selected_contribs)))
                    contrib.generate_amplitudes()