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
import madgraph.iolibs.file_writers as writers

import madgraph.various.misc as misc
from madgraph import MadGraph5Error, InvalidCmd, MG5DIR
import madgraph.core.accessors as accessors
import madgraph.various.banner as banner_mod

from madgraph.iolibs.files import cp, ln, mv

import madgraph.various.misc as misc

_file_path = os.path.split(os.path.dirname(os.path.realpath(__file__)))[0]
logger = logging.getLogger('madgraph.RustExporter')
pjoin = os.path.join


class RustExporterException(MadGraph5Error):
    """ Daughter class specifying an exception occurring within the RustExporter. """
    pass

# ===============================================================================
# RustExporter
# ===============================================================================


class RustExporter(object):
    """Class to take care of exporting Rust source code for generating a low-level
     representation of MadNkLO integrands.
    This class does not inherit from VirtualExporter because it is not responsible
    for exporting one set of amplitudes, but rather a list of MadNkLO integrands."""

    static_template_path = pjoin(_file_path, os.pardir, 'Template', 'rust')
    dynamic_template_path = pjoin(
        _file_path, os.pardir, 'madgraph', 'iolibs', 'template_files', 'rust')

    def __init__(self, cmd_interface, export_options={}):
        """Initialize an Rust exporter with an instance of the main command line interface as well as additional
         export options which provide an access to all necessary information about how the rust backend should
         be exported, together with the list of MadNkLO integrands that must be exported."""

        self.cmd_interface = cmd_interface
        self.export_dir = cmd_interface._export_dir
        self.options = cmd_interface.options
        self.model = cmd_interface._curr_model
        self.export_options = export_options

    def copy_template(self):
        """ Copy the template directories."""

        rust_export_path = pjoin(self.export_dir, 'rust')

        # The overall exporter of MAdNkLO should have already done this copy of the
        # static rust template, so it does not need to be repeated here. But if more static
        # information needs to be added, then it can be done here.
        if not os.path.exists(rust_export_path):
            shutil.copytree(self.static_template_path, self.export_dir)
        # Create an empty integrands directory for hosting integrand implementation here
        # (an empty directory cannot be placed in the static_template_path because it cannot be tracked by git).
        os.makedirs(pjoin(rust_export_path, 'Cards'))
        os.makedirs(pjoin(self.export_dir,'lib','matrix_elements'))

    def export_global_resources(self, all_MEAccessors, repl_dict):
        """ Export global resources independent of any integrand or accessor, basically the base skeleton of the rust
        implementation. """

        # Copy the some static rust template files
        self.copy_template()

        # Save the process output name placeholder
        repl_dict['output_name'] = os.path.basename(self.export_dir)

        # Also set here what is the prefix of the C_bindings matrix element routines
        repl_dict['C_binding_prefix'] = 'C_'

        # Setup the placeholder that will contain the code for instantiating integrands
        repl_dict['instantiate_integrands'] = ''

        # Copy the build.rs with the path to the model file
        rust_writer = writers.RustWriter(
            pjoin(self.export_dir, 'rust', 'madnklo', 'build.rs'), opt='w')
        repl_dict = {}
        repl_dict['model_path'] = 'lib'
        repl_dict['matrix_element_path'] = pjoin('lib', 'matrix_elements')
        rust_writer.writelines(
            open(pjoin(self.dynamic_template_path, 'build.rs')).read(),
            context={}, replace_dictionary=repl_dict
        )

        # Below one should perform the export steps relating to rendering all entries of the
        # all_MEAccessors dictionary exposed to rust as well
        rust_writer = writers.RustWriter(
            pjoin(self.export_dir, 'rust', 'madnklo', 'src', 'matrix_elements.rs'), opt='a')

        # We must then make sure that all accessor are correctly initialised and the underlying resources
        # (such as f2py bindings for example) are compiled. (we may need them to access for instance the list
        # of available
        for (process_key, (me_accessor, pdg_map)) in all_MEAccessors.items():
            all_MEAccessors[process_key] = (me_accessor.__class__.initialize_from_dump(
                                                    me_accessor.generate_dump(), self.export_dir, self.model),pdg_map)
        all_MEAccessors.synchronize()

        for (process_key, (me_accessor, pdg_map)) in all_MEAccessors.items():

            # For we only support exporting matrix element accessors
            if isinstance(me_accessor, accessors.F2PYMEAccessor):
                ME_rs_template = open(pjoin(self.dynamic_template_path,'matrix_element.rs'),'r').read()
                repl_dict = {}
                repl_dict['matrix_element_id']  = me_accessor.get_id()
                repl_dict['matrix_element_lib'] = me_accessor.get_library_name()
                rust_writer.writelines(ME_rs_template, context={}, replace_dictionary=repl_dict)
        
        rust_writer.close()

        return

    def translate_low_level_code(self, low_level_code, function_prefix=''):
        """ Translates the metacode for calling matrix element or computing counterterm into rust code."""

        rust_code = []
        necessary_imports = []
        replacement_dict = {'function_prefix':function_prefix}

        def translate_argument(arg):
            if isinstance(arg, tuple) and arg[0]=='EpsilonExpansion':
                return "epsilonexpansion![%s]"%(
                    ','.join('%d => %s'%(k, translate_argument(v)) for k, v in sorted(arg[1].items(), key=lambda k: k[0]))
                )
            elif isinstance(arg, list):
                return 'Vec![%s]'%(','.join(translate_argument(element) for element in arg))
            else:
                return str(arg)

        def translate_function_call(function_identifier, arguments):
            return '%s(%s)'%(
                function_identifier%replacement_dict,
                ','.join(translate_argument(arg) for arg in arguments)
            )

        for instruction in low_level_code:

            if instruction[0] == 'use':
                replacement_dict[instruction[1]] = instruction[2]
                if instruction[1] == 'ME_library':
                    necessary_imports.append('use crate::MatrixElements::%s;'%instruction[2])
            elif instruction[0] == 'set':
                # Assume immutable for now
                rust_code.append('let %s = %s;'%(instruction[1], translate_argument(instruction[2])))
            elif instruction[0] == 'call':
                rust_code.append('%s;'%translate_function_call(instruction[1],instruction[2]))
            elif instruction[0] == 'call_and_set':
                rust_code.append('let %s = %s;'%(instruction[1], translate_function_call(instruction[2],instruction[3])))
            elif instruction[0] == 'return':
                rust_code.append('%s'%translate_argument(instruction[1]))

        return rust_code, necessary_imports

    def export(self, integrand, repl_dict):
        """ Export one particular integrand. """
        #import pdb
        #pdb.set_trace()

        # The rust exporter currently only support LO integrands
        if integrand.contribution_definition.overall_correction_order.count('N')>0:
            return

        integrand_short_name = '%s_%d' % (
            integrand.get_short_name(), integrand.ID)
        integrand_export_path = pjoin(
            self.export_dir, 'rust', 'madnklo', 'src', 'integrands', integrand_short_name)
        # Create a directory specifically for this integrand
        os.makedirs(integrand_export_path)

        #shutil.copy(
        #    pjoin(self.dynamic_template_path, 'integrand.rs'),
        #    pjoin(integrand_export_path, 'integrand_%s.rs' % (integrand_short_name))
        #)

        # Make sure to first disable accessor caches
        accessors.deactivate_cache()

        # Store here all low-level instructions for performing Matrix Element calls
        ME_calls_instructions = {}
        all_flavor_configurations_per_process = {}

        for i_process, (process_key, (process, mapped_processes)) in enumerate(sorted(integrand.processes_map.items())):

            # First construct the list of
            process_pdgs = process.get_cached_initial_final_pdgs()
            all_processes = [process,]+mapped_processes
            all_flavor_configurations = []
            # The process mirroring is accounted for at the very end only
            for proc in all_processes:
                initial_final_pdgs = proc.get_cached_initial_final_pdgs()
                all_flavor_configurations.append(initial_final_pdgs)
            all_flavor_configurations_per_process[i_process] = all_flavor_configurations

            # The counterterm counter starts at zero for the "physical contribution" (i.e. ME)
            i_CT = 0
            # The i_CT=0 contribution of course then involves only a single call to the ME
            i_call = 0
            assert(not hasattr(integrand, 'counterterms'))
            assert(not hasattr(integrand, 'integrated_counterterms'))

            call_key = (i_process, i_CT, i_call)

            # We can now use the all_flavor_configurations built here to hard-code it in the rust template.
            PS_point = {}
            for i_leg, leg in enumerate(process.get_initial_legs()):
                PS_point[leg.get('number')] = 'p[%d]'%leg.get('number')
            for i_leg, leg in enumerate(process.get_final_legs()):
                PS_point[leg.get('number')] = 'p[%d]'%leg.get('number')
            alpha_s = "alpha_s"
            mu_r = "mu_r"

            # We can now reconstruct what are the necessary ME calls in this sigma call so as to be able to hardcode them
            low_level_calling_code = integrand.all_MEAccessors(
                process, PS_point, alpha_s, mu_r, pdgs=all_flavor_configurations[0], low_level_code_generation=True)

            ME_calls_instructions[call_key] = low_level_calling_code

        # We must now hard-code the collected low-level code into the dynamic rust templates
        # ==================================================================================

        # First output the code for instantiating this particular integrand in the AllIntegrands struct
        integrand_instantiation = []
        integrand_instantiation.append("// Now instantiating integrand '%s'"%integrand_short_name)
        integrand_instantiation.append('all_integrands.all_integrands.insert(%d, Integrand('%integrand.ID)
        instantiation_repl_dict = {}
        instantiation_repl_dict['n_processes'] = len(integrand.processes_map)
        instantiation_repl_dict['n_initial'] = len(integrand.masses[0])
        instantiation_repl_dict['n_final'] = len(integrand.masses[1])
        instantiation_repl_dict['masses'] = "(Vec![%s],Vec![%s])"%(
            ','.join('%s'%m for m in integrand.masses[0]),','.join('%s'%m for m in integrand.masses[1]))
        instantiation_repl_dict['all_flavor_configurations'] = "hashmap![%s]"%(','.join('%d => %s'%
            (i_process, "Vec![%s]"%(','.join('(Vec![%s],Vec![%s])'%(
                ','.join('%d'%pdg for pdg in flavor_config[0]),
                ','.join('%d'%pdg for pdg in flavor_config[1])) for flavor_config in flavor_configs)))
            for i_process, flavor_configs in sorted(all_flavor_configurations_per_process.items(), key=lambda k:k[0])))
        instantiation_repl_dict['run_card'] = "RunCard('Cards/run_card.yaml')"

        integrand_instantiation.append(
            ('%(n_processes)d, %(all_flavor_configurations)s, '+
            '%(n_initial)d, %(n_final)d, %(masses)s, %(run_card)s')%instantiation_repl_dict)
        integrand_instantiation.append(')')
        repl_dict['instantiate_integrands'] += '\n'.join(integrand_instantiation) + '\n'

        # Now we must write out the corresponding integrand evaluator
        integrand_evaluator_repl_dict = {}
        integrand_evaluator_repl_dict['integrand_ID'] = integrand.ID
        ME_calls_lines = []
        evaluator_header = []
        for call_key, low_level_code in sorted(ME_calls_instructions.items(), lambda k: k[0]):
            ME_calls_lines.append('(%d, %d, %d) => {'%call_key)
            rust_code, headers = self.translate_low_level_code(low_level_code, function_prefix=repl_dict['C_binding_prefix'])
            evaluator_header.extend(headers)
            ME_calls_lines.extend(rust_code)
            ME_calls_lines.append('},')
        integrand_evaluator_repl_dict['ME_calls'] = '\n'.join(ME_calls_lines)

        rust_writer = writers.RustWriter(pjoin(integrand_export_path,'integrand_evaluator_%s.rs'%integrand_short_name),opt='w')
        integrand_evaluator_context = {}
        integrand_evaluator_repl_dict['header'] = '\n'.join(evaluator_header)
        rust_writer.writelines(
            open(pjoin(self.dynamic_template_path, 'integrand_evaluator.rs'), 'r').read(),
            context=integrand_evaluator_context, replace_dictionary=integrand_evaluator_repl_dict)
        rust_writer.close()

        # The code below would be an example of how to aggregated all integrand evaluators into a single rust file
        # rust_writer = writers.RustWriter(None)
        # integrand_evaluator_context = {}
        # Headers will be added together in finalize
        # integrand_evaluator_repl_dict['header'] = ''
        # rust_writer.writelines(
        #     open(pjoin(self.dynamic_template_path, 'integrand_evaluator.rs'), 'r').read(),
        #     context=integrand_evaluator_context, replace_dictionary=integrand_evaluator_repl_dict)
        # if 'integrand_evaluators_header' in repl_dict:
        #     repl_dict['integrand_evaluators_header'] += '\n'.join(evaluator_header) + '\n'
        # else:
        #     repl_dict['integrand_evaluators_header'] = '\n'.join(evaluator_header) + '\n'
        # if 'integrand_evaluators_body' in repl_dict:
        #     repl_dict['integrand_evaluators_body'] += rust_writer.pop_content() + '\n'
        # else:
        #     repl_dict['integrand_evaluators_body'] = rust_writer.pop_content() + '\n'

        return

        # ==========================================
        # BELOW IS JUST AN EXAMPLE OF DYNAMIC EXPORT
        # ==========================================

        # Template content:

        # /*
        # * This is the rust implementation of the integrand % (integrand_name)s
        # */
        #
        # ## if (my_conditional_variable > 3) {
        # println!("A");
        # ## } else {
        # println!("B");
        # ## if (my_nested_conditional_boolean) {
        # println!("B1");
        # ## } else {
        # println!("B2");
        # ## }
        # ## }
        #
        # println!("This is inconditional");
        #
        # ## if (any(var=='PrintTHIS' for var in myContextList)) {
        # println!("Some statement dynamically assigned:")
        # % (a_dynamical_statement)
        # s
        # ## }

        # # Create in their the dynamically generate rust source code
        # my_context = {
        #     'my_conditional_variable': 4,
        #     'my_nested_conditional_boolean': True,
        #     'myContextList': ["A sheep", "A tulip", "PrintTHIS"],
        # }
        # my_replacement_dict = {
        #     'integrand_name': integrand_short_name,
        #     'matrix_element_lib': pjoin(self.export_dir, "test"),
        #     'a_dynamical_statement': 'println!("MadNkLO ftw!");'
        # }
        # rust_writer = writers.RustWriter(pjoin(
        #     integrand_export_path, 'integrand_%s.rs' % (integrand_short_name)), opt='w')
        # rust_writer.writelines(
        #     open(pjoin(self.dynamic_template_path, 'integrand.rs')).read(),
        #     context=my_context, replace_dictionary=my_replacement_dict
        # )

    def compile(self):
        """ Compile the rust backend."""
        # TODO Note: this should be done as much as possible using makefile targets so that one can easily recompile
        # the entire distribution manually.
        return

    def finalize(self, all_MEAccessors, all_integrands, repl_dict):
        """Distribute and organize the finalization export of all accessors and integrands. """

        # Now that all low-level rust code has been generated for all integrands and accessor, we can write it on file.

        rust_export_path = pjoin(self.export_dir, 'rust', 'madnklo', 'src')
        # Write the all_integrands instantiator
        rust_writer = writers.RustWriter(pjoin(rust_export_path, 'all_integrands.rs'), opt='w')
        rust_writer.writelines(
            open(pjoin(self.dynamic_template_path, 'integrand_instantiator.rs'), 'r').read(),
            context={}, replace_dictionary=repl_dict)
        rust_writer.close()

        # Also write the yaml equivalent of the run_card
        # Write down there an example of a yaml-formatted run card
        open(pjoin(self.export_dir, 'rust','Cards','run_card.yaml'),'w').write(
            banner_mod.RunCardME7(pjoin(self.export_dir, 'Cards', 'run_card.dat')).export(
                pjoin(self.export_dir, 'Cards', 'run_card.dat'), format='yaml'))