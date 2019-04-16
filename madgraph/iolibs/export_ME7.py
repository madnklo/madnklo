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

"""Methods and classes to export a list of contributions in the ME7_format."""

import glob
import logging
import os
import re
import shutil

import madgraph.core.base_objects as base_objects
import madgraph.core.color_algebra as color
import madgraph.core.helas_objects as helas_objects
import madgraph.core.contributions as contributions
import madgraph.core.subtraction as subtraction
import madgraph.core.accessors as accessors
import madgraph.core.diagram_generation as diagram_generation
import madgraph.interface.common_run_interface as common_run_interface
import aloha as aloha
import models.import_ufo as import_ufo
import models.model_reader as model_reader
import madgraph.iolibs.files as files
import madgraph.iolibs.file_writers as writers
import madgraph.iolibs.template_files as template_files
import madgraph.iolibs.export_v4 as export_v4
import madgraph.iolibs.save_load_object as save_load_object
import madgraph.various.banner as banner_mod
from madgraph import MadGraph5Error, InvalidCmd, MG5DIR
from madgraph.iolibs.files import cp, ln, mv

import madgraph.various.misc as misc

_file_path = os.path.split(os.path.dirname(os.path.realpath(__file__)))[0] + '/'
logger = logging.getLogger('madgraph.ME7Exporter')
pjoin = os.path.join

template_path = pjoin(_file_path,os.pardir,'Template')

#===============================================================================
# ME7Exporter
#===============================================================================
class ME7Exporter(object):
    """Class to take care of exporting a set Contributions in the ME7 format.
    This class does not inherit from VirtualExporter because it is not responsible
    for exporting one set of amplitudes, but rather a list of Contributions."""

    def __init__(self, cmd_interface, noclean, group_subprocesses, export_options={}):
        """Initialize an ME7Exporter with an output path and a list of contributions"""

        self.group_subprocesses = group_subprocesses
        self.cmd_interface      = cmd_interface
        self.contributions      = cmd_interface._curr_contribs
        self.export_dir         = cmd_interface._export_dir 
        self.options            = cmd_interface.options
        self.model              = cmd_interface._curr_model
        self.export_options     = export_options
        # Set some options to default if absent
        if 'ignore_contributions' not in self.export_options:
            self.export_options['ignore_contributions'] = []
        if 'ignore_integrated_counterterms' not in self.export_options:
            self.export_options['ignore_integrated_counterterms'] = []

        # Already initialize the exporter of each contribution
        for contrib in self.contributions:
            contrib.initialize_exporter(
                cmd_interface, noclean, group_subprocesses=group_subprocesses)

        # During the export of a bunch of contributions, store the integrated
        # counterterms that did not find a host contribution, so as to check, *after*
        # all the contributions have been finalized, that the reduced processes of
        # those integrated counterterm is indeed inexistent (i.e. by looking it up in
        # all_ME_accessor dict.)
        self.integrated_counterterms_refused_from_all_contribs = []

    def pass_information_from_cmd(self, cmd_interface):
        """ Pass information from the cmd_interface to this exporter and all contributions."""
        
        # For ME7, apply this to all individual contributions
        for contrib in self.contributions:
            contrib.pass_information_from_cmd(cmd_interface)

    def update_make_opts(self):
        """Synchronize make_opts with the MG5/ME7 options."""

        make_opts = pjoin(self.export_dir, 'Source', 'make_opts')
    
        cpp_compiler = self.options['cpp_compiler'] if self.options['cpp_compiler'] else 'g++'
        fortran_compiler = self.options['fortran_compiler'] if self.options['fortran_compiler'] else 'gfortran'
        f2py_compiler = self.options['f2py_compiler'] if self.options['f2py_compiler'] else 'f2py'

        is_clang = misc.detect_if_cpp_compiler_is_clang(cpp_compiler)
        is_lc    = misc.detect_cpp_std_lib_dependence(cpp_compiler) == '-lc++'

        # list of the variable to set in the make_opts file
        for_update= {'MACFLAG' : '-mmacosx-version-min=10.7' if is_clang and is_lc else '',
                     'STDLIB' : '-lc++' if is_lc else '-lstdc++',
                     'STDLIB_FLAG' : '-stdlib=libc++' if is_lc and is_clang else '',
                     'DEFAULT_CPP_COMPILER' : cpp_compiler,
                     'DEFAULT_F_COMPILER' : fortran_compiler,
                     'DEFAULT_F2PY_COMPILER' : f2py_compiler}
    
        common_run_interface.CommonRunCmd.update_make_opts_full(make_opts, for_update)

    def copy_template(self, model):
        """ Copy the template directories."""

        # The root directory should have already been cleaned if necessary
        if os.path.exists(self.export_dir):
            raise InvalidCmd("The output path '%s' already exists. Clean it first."%self.export_dir)
        
        shutil.copytree(pjoin(template_path,'ME7'),self.export_dir)

        # Generate a .gitignore file so that the process folder is ignored by the versionning system
        gitignore = open(pjoin(self.export_dir,".gitignore"),"w+")
        gitignore.write("*\n")
        gitignore.close()

        # Make sure to have make_opts synchronized
        self.update_make_opts()

        # Forward the request for copying the template to each contribution
        self.contributions.apply_method_to_all_contribs('copy_template', method_args = [model])

    def generate_beam_factorization_contributions_for_correction_order(self, 
                                                 correction_order, n_unresolved_particles):
        """ Generates all beam factorization contributions for the specified correction order
        and number of unresolved particles."""
   
        # Find a representative (VV, RV etc...) contribution of this class whose
        # contribution definition can be used as a template for generating the beam
        # factorization ones, Skip if None is found.
        candidates = [c for c in self.contributions.get_contributions_of_order(
            'N'*correction_order+'LO') if c.contribution_definition.n_unresolved_particles==
            n_unresolved_particles]

        if len(candidates)==0:
            logger.critical('No host contribution was found to be template of the beam factorization'+
                ' contribution of peturbative order %s with %d unresolved final state%s.\n'%(
                    'N'*correction_order+'LO', n_unresolved_particles, 's' if n_unresolved_particles > 1 else '')+
                'Results from this output are therefore incomplete and unphysical.')
            return []
        
        # Now split the candidates according to their process_IDs
        candidates_per_process_IDs = {}
        for contrib in candidates:
            proc_ID = contrib.contribution_definition.process_definition.get('id')
            if proc_ID in candidates_per_process_IDs:
                candidates_per_process_IDs[proc_ID].append(contrib)
            else:
                candidates_per_process_IDs[proc_ID] = [contrib,]
        
        # Now generate the beam factorization contributions independently for each group
        # of contributions with one particular process ID.
        beam_factorization_contributions_added = []
        for process_ID, template_contributions in candidates_per_process_IDs.items():
            beam_factorization_contributions_added.extend(
                self.generate_beam_factorization_contributions_for_correction_order_and_process_ID(
                    correction_order, n_unresolved_particles, process_ID, template_contributions[0]))

        return beam_factorization_contributions_added

    def generate_beam_factorization_contributions_for_correction_order_and_process_ID(self, 
                    correction_order, n_unresolved_particles, process_ID, template_contribution):
        """ Generates all beam factorization contributions for the specified correction order
        and number of unresolved particles as well as given process_ID. The template contribution
        given in argument is the non-beam-factorized contribution of type V, RV, VV, etc...
        whose contribution definition should serve as a template for generating the beam
        factorization contributions."""

        beam_types = template_contribution.contribution_definition.get_beam_types()

        # Deduce what max_correction order can be put in beam factorization terms here
        beam_factorization_max_order = correction_order-n_unresolved_particles

        # Now generate the basic options to pass to the constructor of the beam factorization
        # contribution
        contrib_def_options = {
            'overall_correction_order'    : template_contribution.contribution_definition.overall_correction_order,
            'correction_order'            : template_contribution.contribution_definition.correction_order,
            'correction_couplings'        : template_contribution.contribution_definition.correction_couplings,
            'squared_orders_constraints'  : template_contribution.contribution_definition.squared_orders_constraints,
            'beam_types'                  : beam_types,
            'n_unresolved_particles'      : n_unresolved_particles,
            'beam_factorization_active'   : (False, False),
            'correlated_beam_convolution' : False
        }

        # When constructing the beam factorization contributions, we will assign 
        # the only argument 'process_definition' to an empty one specifying only limited
        # features like the process ID since in this case the processes corresponding 
        # to these beam factorization contributions will effectively be pulled from 
        # the various other existing contributions. 
        # Most of the process definition attributes are anyway only useful at generation 
        # time.
        beam_factorization_contributions = []
        if not beam_types[0] is None:
            contrib_def_options['beam_factorization_active'] = (True, False)
            contrib_def_options['n_loops']                   = correction_order-n_unresolved_particles-1
            dummy_process_definition = base_objects.ProcessDefinition({
                'model'     : template_contribution.contribution_definition.process_definition.get('model'),
                'id'        : process_ID,
                'n_loops'   : correction_order-n_unresolved_particles-1
            })
            beam_factorization_contributions.append( contributions.Contribution(
                base_objects.ContributionDefinition(dummy_process_definition,
                                              **contrib_def_options ),self.cmd_interface) )

        if not beam_types[1] is None:
            contrib_def_options['beam_factorization_active'] = (False, True)
            contrib_def_options['n_loops']                   = correction_order-n_unresolved_particles-1
            dummy_process_definition = base_objects.ProcessDefinition({
                'model'     : template_contribution.contribution_definition.process_definition.get('model'),
                'id'        : process_ID,
                'n_loops'   : correction_order-n_unresolved_particles-1
            })
            beam_factorization_contributions.append( contributions.Contribution(
                base_objects.ContributionDefinition(dummy_process_definition,
                                              **contrib_def_options ),self.cmd_interface) )
        if (not beam_types[0] is None) and (not beam_types[1] is None) and \
                                                         beam_factorization_max_order >= 2:
            contrib_def_options['beam_factorization_active'] = (True, True)
            contrib_def_options['n_loops']                   = correction_order-n_unresolved_particles-2
            dummy_process_definition = base_objects.ProcessDefinition({
                'model'     : template_contribution.contribution_definition.process_definition.get('model'),
                'id'        : process_ID,
                'n_loops'   : correction_order-n_unresolved_particles-2
            })
            beam_factorization_contributions.append( contributions.Contribution(
                base_objects.ContributionDefinition(dummy_process_definition,
                                              **contrib_def_options ),self.cmd_interface) )
        
        # Now add the soft beam factorization contributions necessary for absorbing the colorful
        # soft integrated CT when using a mapping that recoils equally against both initial states.
        if ((not beam_types[0] is None) and (not beam_types[1] is None)) and \
           self.options['subtraction_currents_scheme'] in subtraction._currents_schemes_requiring_soft_beam_factorization:
            contrib_def_options['beam_factorization_active']   = (True, True)
            contrib_def_options['n_loops']                     = correction_order-n_unresolved_particles-1
            contrib_def_options['correlated_beam_convolution'] = True            
            dummy_process_definition = base_objects.ProcessDefinition({
                'model'     : template_contribution.contribution_definition.process_definition.get('model'),
                'id'        : process_ID,
                'n_loops'   : correction_order-n_unresolved_particles-1
            })
            beam_factorization_contributions.append( contributions.Contribution(
                base_objects.ContributionDefinition(dummy_process_definition,
                                              **contrib_def_options ),self.cmd_interface) )

        if len(beam_factorization_contributions)==0:
            return []

        # Now we fish out all lower order contributions and their processes and add them
        # to the contributions we instantiated above:
        for lower_order in range(correction_order):
            for contrib in self.contributions.get_contributions_of_order('N'*lower_order+'LO'):
                # Make sure it matches the number of unresolved particles
                if contrib.contribution_definition.n_unresolved_particles != n_unresolved_particles:
                    continue
                # And make sure it has the correct process ID
                if contrib.contribution_definition.process_definition.get('id') != process_ID:
                    continue
                # Now deduce what max_correction order can be put in beam factorization terms 
                # for this particular contribution
                beam_factorization_order = correction_order - lower_order

                # Now add the processes from contrib to the beam factorization contributions
                # that we just instantiated.     
                for beam_factorization_contrib in beam_factorization_contributions:
                    beam_factorization_contrib.add_beam_factorization_processes_from_contribution(
                                                         contrib, beam_factorization_order)

        # Return all contributions instantiated here
        return beam_factorization_contributions

    def generate_all_beam_factorization_contributions(self):
        """ Generates all necessary additional beam factorization contributions, also
        assigning them them a reference to the already generated attributes of the original
        contributions."""
        
        # If no contribution has an "active" beam then simply return here immediately
        # Actually, nothing wrong would happen if we were to continue nonetheless, but 
        # some info/warning message could be sent which would be awkward for the user to
        # understand if, say, he is doing some e+ e- process. We therefore opt for returning
        # here immediately.
        if all( c.contribution_definition.beam_factorization==
                         {'beam_one': None, 'beam_two': None} for c in self.contributions):
            return

        max_correction_order = self.contributions.get_max_correction_order().count('N')
        
        beam_factorization_contributions = []
        # We will now add the beam factorization contributions at each perturbative order
        for correction_order in range(max_correction_order):
            for n_unresolved_particles in range(correction_order+1):
                # For each possible number of unresolved particles, create a maximum of 
                # three new beam_factorization contribution:
                # At NLO, we would create the contributions VF1 and VF2. 
                # At NNLO RV would be use as a template for RVF1 and RVF2 and VV for VVF1, VVF2 and VVF1F2 
                # etc...   
                # In addition, when using colorful with a soft mapping that recoils against
                # initial state, we will generate an extra contribution for each correction order
                # and number of unresolved particles of the form:
                # BS at NLO
                # VS, RS at NNLO, etc...
                # where both beam are convoluted with an identical fraction xi.
                # ------------------------------------------------------------------------------
                # /!\ Notice that depending on the specifics of the implementation of colorful
                # with ISR at NNLO, it may be that contributions such as VVF1S may be needed.
                # These will not be generated for now, but it can be accommodated in the future
                # by generating the soft beam convolution contributions *after* and *on top of*
                # the beam factorization contributions already generated.
                # ------------------------------------------------------------------------------
                beam_factorization_contributions.extend(
                    self.generate_beam_factorization_contributions_for_correction_order(
                                               correction_order+1, n_unresolved_particles))

        # Add the new contributions generated unless explicitly asked to be ignored by the user
        for bf_contribution in beam_factorization_contributions:
            if bf_contribution.short_name() not in self.export_options['ignore_contributions']:
                self.contributions.append(bf_contribution)
            else:
                logger.warning("User explicitly asked to remove contribution "+
                    "'%s'. This can potentially yield incorrect results."%bf_contribution.short_name())
        self.contributions.sort_contributions()
    
    def export(self, nojpeg, args=[]):
        """ Distribute and organize the export of all contributions. """

        # Forward the export request to each contribution
        export_return_values = self.contributions.apply_method_to_all_contribs('export', 
            method_args = [],
            method_opts = {
                'nojpeg':nojpeg, 
                'group_processes':self.group_subprocesses,
                'args':args}
        )

        # Now generate all necessary additional beam_factorization contributions
        logger.info('Adding beam factorization contributions...')
        self.generate_all_beam_factorization_contributions()

        # Now, we can perform subtraction

        # Pass the options to specify for each contrib as a callable that must be evaluated
        # for each contrib
        def method_options(contrib):
            shortname = contrib.short_name()
            to_ignore = self.export_options['ignore_integrated_counterterms']
            return {'nojpeg':nojpeg, 
                    'group_processes':self.group_subprocesses,
                    'ignore_integrated_counterterms': ('all' in to_ignore) or (shortname in to_ignore),
                    'args':args}

        return_values = self.contributions.apply_method_to_all_contribs('subtract', 
            method_args = [],
            method_opts = method_options)

        # Now gather all integrated counterterms generated during the export and dispatch
        # them to the right contributions (i.e. the integrated counterterms from the 
        # real-emission contribution belong to the virtual contribution.)
        for contribution, return_dict in return_values:
            if 'integrated_counterterms' in return_dict:
                self.distribute_integrated_counterterms(contribution, 
                                                    return_dict['integrated_counterterms'])
        
    def distribute_integrated_counterterms(self, contribution_origin, integrated_counterterms):
        """Analyses the integrated counterterms from the list provided in argument and 
        coming from the specified contribution (presumably some real-emission type of
        contributions) and assign them to the correct contribution
        (typically of virtual origin)."""
                
        # Gather which contribution will receive the integrated counterterm for a given
        # number of loops and unresolved legs, and of course a given process_defining_ID
        # so that contributions from different 'add process' commands don't get mangled.
        # Note that in principle one could do a look-up of the processes in the process_map
        # of each contributions, but this is unnecessarily slow and complicated; simply
        # assigning the routing map of the distribution of the integrated counterterms
        # based on the (proc_ID, n_loops, n_unresolved, beam_one_active, beam_two_active, 
        #               correlated_beam_convolution) 
        # is enough.
        # However, when subprocesses grouping is active (especially mirroring), the initial
        # states can be swapped, such that an integrated counterterm naively requiring
        # beam one to be convoluted must actually be placed in BF2. Example:
        #
        # Real process 'g d~ > z d~' yields '[C(1,3)]' with reduced process 'd d~ > z' 
        # which must indeed be placed in the contribution BF1 'd d~ > z'
        # However, Real process 'g d > z d' also yields '[C(1,3)]' but with the reduced
        # process 'd~ d > z' which must this time be placed in the contribution BF*2* 'd d~ > z'
        #
        # For this reason, when grouping subprocesses, we will attempt to place the integrated
        # CT in both BF1 and BF2 and the implementation of 'add_integrated_counterterm' is such
        # that only the right beam-factorization contribution will accept it.
        def key_string(key):
            """ Return a nice string of the routing key in this function."""
            if self.group_subprocesses:
                return ("(proc_ID=%d, n_loops=%d, n_unresolved=%d, n_beams_active=%s, "+
                        "correlated_beam_convolution=%s)")%key
            else:
                return ("(proc_ID=%d, n_loops=%d, n_unresolved=%d, beam_one_active=%s, "+
                        "beam_two_active=%s, correlated_beam_convolution=%s)")%key
            
        routing_map = {}
        for contribution in self.contributions:
            n_active_beams = 0
            if contribution.contribution_definition.is_beam_active('beam_one'):
                n_active_beams += 1
            if contribution.contribution_definition.is_beam_active('beam_two'):
                n_active_beams += 1
            if not self.group_subprocesses:
                key = (
                    contribution.contribution_definition.process_definition.get('id'),
                    contribution.contribution_definition.n_loops,
                    contribution.contribution_definition.n_unresolved_particles,
                    contribution.contribution_definition.is_beam_active('beam_one'),
                    contribution.contribution_definition.is_beam_active('beam_two'),
                    contribution.contribution_definition.correlated_beam_convolution
                )
            else:
                key = (
                    contribution.contribution_definition.process_definition.get('id'),
                    contribution.contribution_definition.n_loops,
                    contribution.contribution_definition.n_unresolved_particles,
                    n_active_beams,
                    contribution.contribution_definition.correlated_beam_convolution
                )
            if key in routing_map:
                # Only trigger the warning if not grouping subprocesses, otherwise it is
                # expected that BF1 and BF2 share the same key, by design
                if not self.group_subprocesses or n_active_beams!=1:
                    logger.warning("The two contributions:\n    %s\nand\n    %s\n"%(
                        routing_map[key][0].nice_string(),contribution.nice_string())+
                        ("share the same  key %s.\n")%key_string(key)+
                        "All integrated counterterms will be placed in the first of the matching"+
                        " contributions that accepts the integrated counterterm.")
                routing_map[key].append(contribution)
            else:
                routing_map[key] = [contribution,]

        # Only warn once about issues occurring during the dispatching of integrated
        # counterterms
        warned = False
        for integrated_CT_properties in integrated_counterterms:
            # Extract quantities from integrated_counterterm_properties
            counterterm = integrated_CT_properties['integrated_counterterm']
            proc_def_ID = counterterm.process.get('id')
            
            # Decide whether this counterterm need to be hosted in a contribution with
            # beam factorization
            necessary_beam_convolutions = counterterm.get_necessary_beam_convolutions()
            beam_one_convolution = 'beam_one' in necessary_beam_convolutions
            beam_two_convolution = 'beam_two' in necessary_beam_convolutions
            n_active_beams = 0
            if beam_one_convolution:
                n_active_beams += 1
            if beam_two_convolution:
                n_active_beams += 1

            # The integrated counterterm belongs to contribution whose defining number of
            # loops is the sum of the number of loops in the reduced process and the 
            # integrated current *plus* the total number of unresolved legs of that structure.
            # For instance, 
            #   > the integrated counterterms of RV belong in VV
            #   > the integrated *single-unresolved* counterterms of RR belong in RV, etc...
            #
            # And beam factorization terms like BornME * F^{(0)}_1 * F^{(0)}_2 are hosted
            # in contributions with n_unresolved = 0 and n_loops = 2 (since they are akin 
            # to VV)
            n_loops = counterterm.n_loops_in_host_contribution()

            # Then the number of unresolved particle of the contribution that receives
            # this counterterm should be the number of unresolved emission of the
            # originating contributions minus the number of unresolved legs in the integrated
            # current.
            n_unresolved = contribution_origin.contribution_definition.n_unresolved_particles - \
                                                             counterterm.count_unresolved()

            # Checks whether this integrated counterterm requires a host contribution featuring
            # correlated convolution of the beams (i.e. 'BS', 'VS', etc... contributions).
            # This is the case only for integrated counterterms with pure-soft *BeamCurrent* that originated
            # from a colorful subtraction currents scheme with mappings such as ppToOneWalker that
            # recoils the soft momentum equally against both initial-state beams.
            correlated_beam_convolution = counterterm.does_require_correlated_beam_convolution()

            if not self.group_subprocesses:
                key = ( proc_def_ID, n_loops, n_unresolved, 
                        beam_one_convolution, beam_two_convolution, correlated_beam_convolution)
            else:
                # As explained before, we must always engineer BF1 and BF2 (similarly at higher orders)
                # as potential matches when grouping subprocesses (with mirroring), so in this case
                # we build the routing key only based off the number of active beam factorization.
                key = ( proc_def_ID, n_loops, n_unresolved, n_active_beams, correlated_beam_convolution)
                
            # This missing contribution can happen if for example the user explicitly 
            # disabled some contributions, typically the virtual.
            # For now we simply skip the incorporation of this integrated counterterms and
            # issue a critical warning, but eventually one can think of creating an 
            # ad-hoc "dummy" contribution to contain those integrated counterterms.
            if key not in routing_map:
                msg = ("Could not find a contribution with key '%s'"%key_string(key)+" to host the"+
                    " integrated counterterm:\n%s\n with key:\n%s"%(counterterm.nice_string(),key_string(key))+
                    "\nIt will therefore be skipped making the ensuing results unphysical and wrong.\n"+
                    "Available keys are:\n%s"%('\n'.join(key_string(k) for k in routing_map)))
                if __debug__:
                    if not warned:
                        logger.critical(msg)
                        warned = True
                        logger.critical("Further occurrences of this warning will now be suppressed.")
                    continue
                else:
                    raise MadGraph5Error(msg)

            found_contribution_host = False
            for contribution_candidate in routing_map[key]:
                # Try to add this integrated counterterm to all candidate contributions
                # and stop at the first one that accepts it based on the processes it
                # contains. Complain if no contribution accepted the counterterm.
                if contribution_candidate.add_integrated_counterterm(integrated_CT_properties):
                    found_contribution_host = True
                    break
            if not found_contribution_host:
                # add this contribution with no host to the list and after all contributions
                # will be finalized, this exporter will check that their reduced processes
                # were indeed inexistent.
                self.integrated_counterterms_refused_from_all_contribs.append(integrated_CT_properties)
                continue
        
        # Finally, the initial state collinear counterterms may be combined because the
        # BeamCurrents that implements them return the full backward evolution flavor
        # matrix, so that several initial state collinear counterterms that differ only
        # by the resulting backward evolved flavor (e.g. 'd g <- d' and 'g d~ <- d') could
        # be combined. This option is not implemented as of now, so the function below will
        # do nothing.
        for contribution in self.contributions:
            contribution.combine_initial_state_counterterms()
                
    def copy_model_resources(self):
        """Make the copy/symbolic links"""
        model_path = pjoin(self.export_dir, 'Source','MODEL')
        if os.path.exists(pjoin(model_path, 'ident_card.dat')):
            mv(pjoin(model_path,'ident_card.dat'), pjoin(self.export_dir, 'Cards'))
        cp(pjoin(model_path, 'param_card.dat'), pjoin(self.export_dir, 'Cards'))
        mv(pjoin(model_path, 'param_card.dat'), pjoin(self.export_dir, 'Cards','param_card_default.dat'))


    def dump_ME7(self, all_MEAccessors, all_integrands):
        """ Dumps all necessary information in order to bootstrap the ME7 interface.
        It is mostly all contained in all_MEAccessors and all_integrands."""
        
        # Simply retrieve this information from one process of the first contribution.
        n_initial = self.contributions[0].get_processes_map().values()[0][0].get_ninitial()
        
        # For now the db is a rude pickle file, but we might improve this to an actual DB eventually
        integrand_dumps = [integrand.generate_dump() for integrand in all_integrands]
        import madgraph
        for itg_dump in integrand_dumps:
            itg_dump['class'] = 'madgraph.integrator.ME7Integrands.ME7_integrands.ME7CythonIntegrand_B'
        ME7_dump = {
            'all_MEAccessors' : all_MEAccessors.generate_dump(),
            'all_integrands'  : [integrand.generate_dump() for integrand in all_integrands],
#            'all_integrands'  : integrand_dumps,
            'model_name'      : 'ME7_UFO_model_%s'%self.model.get('name'),
            'model_with_CMS'  : self.options['complex_mass_scheme'],
            'n_initial'       : n_initial,
            'subtraction_currents_scheme': self.options['subtraction_currents_scheme']
        }
        
        save_load_object.save_to_file(pjoin(self.export_dir,'MadEvent7.db'), ME7_dump)

    def create_run_card(self):
        """ Create the run card."""
        
        run_card = banner_mod.RunCardME7()

        history = ''
        processes = [[v[0] for v in contrib.get_processes_map().values()] for contrib in self.contributions]
        proc_characteristic = {
            'ninitial':processes[0][0].get_ninitial(), 
            'loop_induced': len(self.contributions.get_loop_induced_contributions()), 
            'colored_pdgs': range(1,7)+[21]}

        run_card.create_default_for_process(proc_characteristic, history, processes)

        run_card.write(pjoin(self.export_dir, 'Cards', 'run_card.dat'), 
            template=pjoin(self.export_dir, 'Cards', 'run_card.dat'), python_template=True )
        run_card.write(pjoin(self.export_dir, 'Cards', 'run_card_default.dat'), 
            template=pjoin(self.export_dir, 'Cards', 'run_card.dat'), python_template=True )

    def compile(self):
        """ Compile all contributions and the global ME7 resources (e.g. MODEL)"""
        
        # Compile the MODEL first
        if os.path.isdir(pjoin(self.export_dir,'Source','MODEL')):
            logger.info("Compiling global ME7 Model")
            misc.compile(arg=['../../lib/libmodel.a'], cwd=pjoin(self.export_dir,'Source','MODEL'), mode='fortran')
        
        # Compile all contributions
        self.contributions.apply_method_to_all_contribs('compile', log='Compiling')

    def finalize(self, flaglist, interface_history):
        """Distribute and organize the finalization of all contributions. """
        
        # Make sure contributions are sorted at this stage
        # It is important to act on LO contributions first, then NLO, then etc...
        # because ME and currents must be added to the ME_accessor in order since there
        # are look-up operations on it in-between
        self.contributions.sort_contributions()

        # Save all the global couplings to write out afterwards
        global_wanted_couplings = []
        # Forward the finalize request to each contribution
        for contrib in self.contributions:
            # Must clean the aloha Kernel before each aloha export for each contribution
            aloha.aloha_lib.KERNEL.clean()
            wanted_couplings_to_add_to_global = contrib.finalize(
                flaglist=flaglist, interface_history=interface_history)
            global_wanted_couplings.extend(wanted_couplings_to_add_to_global)

        # Generate the global ME7 MODEL
        if global_wanted_couplings:
            output_dir=pjoin(self.export_dir, 'Source', 'MODEL')
            # Writing out the model common to all the contributions that can share it
            model_export_options = {
                'complex_mass'  : self.options['complex_mass_scheme'],
                'export_format' : 'madloop', # So as to have access to lha_read_mp.f
                'mp'            : True,
                'loop_induced'  : False }
            model_builder = export_v4.UFO_model_to_mg4(self.model, output_dir, model_export_options)
            model_builder.build(global_wanted_couplings)
        
        # Now possibly add content to the pool of global ME7 resources before removing superfluous files
        # and linking to necessary global ME7 resources
        for contrib in self.contributions:
            contrib.add_content_to_global_ME7_resources(self.export_dir)
            contrib.remove_superfluous_content()
            contrib.link_global_ME7_resources(self.export_dir)

        # Create the run_card
        self.create_run_card()

        # Add the cards generated in MODEL to the Cards directory
        self.copy_model_resources()
        # Now link the Sources files within each contribution
        for contrib in self.contributions:
            contrib.make_model_symbolic_link()
        
        # Copy the UFO model to the global ME7 resources Source directory
        ME7_ufo_path = pjoin(self.export_dir,'Source','ME7_UFO_model_%s'%os.path.basename(self.model.get('modelpath')))
        shutil.copytree(self.model.get('modelpath'), ME7_ufo_path)
        # And clear compiled files in it
        for path in misc.glob(pjoin(ME7_ufo_path,'*.pkl'))+misc.glob(pjoin(ME7_ufo_path,'*.pyc')):
            os.remove(path)
        
        # Now generate all the ME accessors and integrand.
        # Notice that some of the information provided here (RunCard, ModelReader, root_path, etc...)
        # can and will be overwritten by the actualized values when the ME7Interface will be launched.
        # We provide it here just so as to be complete.

        # Obtain all the Accessors to the Matrix Element and currents made available in this process output
        all_MEAccessors = accessors.MEAccessorDict()
        for contrib in self.contributions:
            contrib.add_ME_accessors(all_MEAccessors, self.export_dir)
        
        # Now make sure that the integrated counterterms without any contribution host
        # indeed have a non-existent reduced process.
        contributions.Contribution_V.remove_counterterms_with_no_reduced_process(
                   all_MEAccessors, self.integrated_counterterms_refused_from_all_contribs)

        # Check there is none left over after this filtering
        if len(self.integrated_counterterms_refused_from_all_contribs)>0:
            counterterm_list = (
                ct['integrated_counterterm'].nice_string()
                for ct in self.integrated_counterterms_refused_from_all_contribs )
            # These integrated counterterms should in principle been added
            msg = "The following list of integrated counterterm are in principle non-zero"
            msg += " but could not be included in any contributions generated:\n"
            msg += '\n'.join(counterterm_list)
            msg += "\nResults generated from that point on are likely to be physically wrong."
            if __debug__:
                logger.critical(msg)
            else:
                raise MadGraph5Error(msg)

        # Now generate all the integrands from the contributions exported
        all_integrands = []
        run_card = banner_mod.RunCardME7(pjoin(self.export_dir,'Cards','run_card.dat'))

        # We might want to recover whether prefix was used when importing the model and whether
        # the MG5 name conventions was used. But this is a detail that can easily be fixed later.
        modelReader_instance = import_ufo.import_model(
            pjoin(self.export_dir,'Source','ME7_UFO_model_')+self.model.get('name'),
            prefix=True,
            complex_mass_scheme=self.options['complex_mass_scheme'] )
        modelReader_instance.pass_particles_name_in_mg_default()
        modelReader_instance.set_parameters_and_couplings(
                param_card = pjoin(self.export_dir,'Cards','param_card.dat'), 
                scale=run_card['scale'], 
                complex_mass_scheme=self.options['complex_mass_scheme'])

        for contrib in self.contributions:
            all_integrands.extend(
                contrib.get_integrands(
                    modelReader_instance, run_card, all_MEAccessors, self.options ) )
   
        # And finally dump ME7 output information so that all relevant objects
        # can be reconstructed for a future launch with ME7Interface.
        # Normally all the relevant information should simply be encoded in only:
        #  'all_MEAccessors' and 'all_integrands'.
        self.dump_ME7(all_MEAccessors, all_integrands)
        
        # Finally, for future convenience it may sometimes be desirable to already compile 
        # all contributions and global ME7 resources (e.g. MODEL) as followed.
        # By default however, we don't do that and this will instead be done at the launch time.
        #logger.info('Compilation of the process output.')
        #logger.info('It can be interrupted at any time,'+
        #                 ' in which case it would be automatically resumed when launched.')
        #self.compile()

        return
        ###################################################################################################
        ###
        ###  WARNING THE CODE BELOW IS JUST FOR TESTING PURPOSES AND CORRESPONDS TO RUNNING THE INTEGRATION
        ###  RIGHT AWAY AND NOT WITHIN THE ME7 INTERFACE>
        ###
        ###################################################################################################
        
        import madgraph.interface.ME7_interface as ME7_interface
        # Test the reconstruction of the ME7 output instances
        ME7_dump = save_load_object.load_from_file(pjoin(self.export_dir,'MadEvent7.db'))
        all_MEAccessors = ME7_dump['all_MEAccessors']['class'].initialize_from_dump(
                                                ME7_dump['all_MEAccessors'], root_path = self.export_dir)
        all_integrands = [integrand_dump['class'].initialize_from_dump(integrand_dump,
                       modelReader_instance, run_card, all_MEAccessors, self.options
                                                    ) for integrand_dump in ME7_dump['all_integrands']]
        model_name     = ME7_dump['model_name']
        model_with_CMS = ME7_dump['model_with_CMS']
        
        # This is now just for gigs. Integrate that beast!
        # Of course, what should really happen is that the users starts a ME7_interface, that 
        # bootstraps from the dump above and starts the integration below with lunch.
        # So this is really just for testing purposes.
        import madgraph.integrator.integrators as integrators
        integrator_naive = integrators.SimpleMonteCarloIntegrator(all_integrands,
            **{'n_iterations'            : 10,
               'n_points_per_iterations' : 100,
               'accuracy_target'         : None,
               'verbosity'               : 1 }
            )
        import madgraph.integrator.pyCubaIntegrator as pyCubaIntegrator        
        integrator_vegas = pyCubaIntegrator.pyCubaIntegrator(all_integrands, 
            **{'algorithm' : 'Vegas', 
               'verbosity' : 1,
               'seed'      : 3,
               'target_accuracy' : 1.0e-3,
               'n_start'   : 1000,
               'n_increase': 500,
               'n_batch'   : 1000,
               'max_eval'  : 100000,
               'min_eval'  : 0})

        # Now run them all!
        for integrator in [integrator_naive, integrator_vegas]:
            xsec, error = integrator.integrate()
            logger.info("="*100)
            logger.info('{:^100}'.format("\033[92mCross-section for process output '%s' with integrator '%s':\033[0m"
                                                                 %(self.export_dir, integrator.get_name())))
            logger.info('{:^100}'.format("\033[94m%.5e +/- %.2e [pb]\033[0m"%(xsec, error)))
            logger.info("="*100+"\n")
