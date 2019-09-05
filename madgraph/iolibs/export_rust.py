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

import logging
import os
import re
import shutil
import subprocess


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

    def __init__(self, cmd_interface, subtraction_module, export_options={}):
        """Initialize an Rust exporter with an instance of the main command line interface as well as additional
         export options which provide an access to all necessary information about how the rust backend should
         be exported, together with the list of MadNkLO integrands that must be exported."""

        self.cmd_interface = cmd_interface
        self.export_dir = cmd_interface._export_dir
        self.options = cmd_interface.options
        self.model = cmd_interface._curr_model
        self.subtraction_module = subtraction_module
        self.export_options = export_options

        # This dictionary will contain the various information for each library to be built and linked
        self.build_info = {
            'helas' : [], # entries of list are '{'lib_name' : ... , 'source_dir' : ... ,
                          #                       'makefile_env_variables' : {}, 'makefile_target' : }
            'matrix_elements' : [] # same as above,
        }

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
        os.makedirs(pjoin(rust_export_path, 'madnklo', 'src', 'integrands'))
        os.makedirs(pjoin(self.export_dir,'lib','matrix_elements'))
        os.makedirs(pjoin(self.export_dir,'lib','helas'))

    def write_param_card_rust_implementation(self, rust_writer, model):
        """ Write a param_card.rs rust implementation that holds all particle masses and width, as well as alpha_s and mu_r."""

        replacement_dict = {}
        all_parameter_lines = []
        for param_name in model.get('parameter_dict'):
            all_parameter_lines.append('\tpub %s : %s,'%(param_name, 'f64'))

        replacement_dict['model_parameters'] = '\n'.join(all_parameter_lines)
        rust_writer.writelines(
            open(pjoin(self.dynamic_template_path, 'param_card.rs'), 'r').read(),
            context={}, replace_dictionary=replacement_dict)

    def export_global_resources(self, all_contributions, all_MEAccessors, repl_dict):
        """ Export global resources independent of any integrand or accessor, basically the base skeleton of the rust
        implementation. The big list of all_contributions is mostly passed so that we can list the necessary
        library libdhelas to be linked."""

        # Copy the some static rust template files
        self.copy_template()

        # Fill in the information necessary for building the makefile of rust resources
        for contrib in all_contributions:
            export_base_dir = os.path.basename(contrib.export_dir)
            # Contributions like BF, BS and similar do no correspond to a ME output directory
            # and therefore do no need their helas library to be output
            if export_base_dir in [None,'None']:
                continue
            self.build_info['helas'].append({
                'lib_name' : 'helas_%s'%export_base_dir,
                'source_dir' : '$(PROC_ROOT)/%s/Source/DHELAS'%export_base_dir,
                'makefile_env_variables' : {
                    'LIBDIR':'../../../lib/helas/',
                    'LIBRARY': 'lib%s.a'%export_base_dir},
                'makefile_target' : '../../../lib/helas/lib%s.a'%export_base_dir
            })

        # Save the process output name placeholder
        repl_dict['output_name'] = os.path.basename(self.export_dir)

        # Also set here what is the prefix of the C_bindings matrix element routines
        repl_dict['C_binding_prefix'] = 'c_'

        # Setup the placeholder that will contain the code for instantiating integrands
        repl_dict['instantiate_integrands'] = ''
        repl_dict['instantiate_integrands_header'] = ''
        repl_dict['integrands_include'] = ''

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

        # Write the model card param_card.rs
        self.write_param_card_rust_implementation(
            writers.RustWriter(pjoin(self.export_dir, 'rust', 'madnklo', 'src', 'param_card.rs'), opt='w'),
            self.model
        )

        # Below one should perform the export steps relating to rendering all entries of the
        # all_MEAccessors dictionary exposed to rust as well
        rust_writer = writers.RustWriter(
            pjoin(self.export_dir, 'rust', 'madnklo', 'src', 'matrix_elements.rs'), opt='a')

        # First add the HELAS linking
        helas_linking_header = []
        for contrib in all_contributions:
            # Ignore contributiosn with no HELAS directories like BFi or BS ones.
            if contrib.export_dir in [None, 'None']:
                continue
            helas_linking_header.append('\n'.join([
                '#[link(name = "%s", kind = "static")]'%os.path.basename(contrib.export_dir),
                'extern "C" {}'
            ]))
        rust_writer.write('\n%s\n'%('\n'.join(helas_linking_header)))

        # We must then make sure that all accessor are correctly initialised and the underlying resources
        # (such as f2py bindings for example) are compiled. (we may need them to access for instance the list
        # of available
        for (process_key, (me_accessor, pdg_map)) in all_MEAccessors.items():
            all_MEAccessors[process_key] = (me_accessor.__class__.initialize_from_dump(
                                                    me_accessor.generate_dump(), self.export_dir, self.model),pdg_map)
        all_MEAccessors.synchronize()

        exported_matrix_element_ids = set()
        subtraction_current_accessor_id = 0
        for (process_key, (me_accessor, pdg_map)) in all_MEAccessors.items():

            # For the rust low-level output each subtraction current accessor needs a unique ID.
            # This unique ID is automatically generated for Matrix element accessors, but we need to issue it here
            # for subtraction accessor.
            if isinstance(me_accessor, accessors.SubtractionCurrentAccessor):
                subtraction_current_accessor_id += 1
                me_accessor.id = subtraction_current_accessor_id
            if me_accessor.get_id() in exported_matrix_element_ids:
                continue
            exported_matrix_element_ids.add(me_accessor.get_id())

            # Logic for rust export of low level matrix element interfaces
            lib_names_to_be_exported = set([])
            if isinstance(me_accessor, accessors.F2PYMEAccessor):
                ME_rs_template = open(pjoin(self.dynamic_template_path,'matrix_element.rs'),'r').read()
                repl_dict = {}
                repl_dict['matrix_element_id']  = me_accessor.get_id()
                repl_dict['matrix_element_lib'] = me_accessor.get_library_name()
                repl_dict['prefix'] = open( pjoin(self.export_dir, '%s/SubProcesses/P%s/proc_prefix.txt'%(
                                                    me_accessor.proc_dir, me_accessor.proc_name)),'r').read().lower()
                if not (me_accessor.get_library_name() in lib_names_to_be_exported):
                    lib_names_to_be_exported.add(me_accessor.get_library_name())

                    # The relative path fo the loop directory is different for the tree-level and MadLoop accessors
                    # Also all SubProcesses of directory are combined within one single directory in the case of
                    # MadLoop output.
                    if isinstance(me_accessor, accessors.F2PYMEAccessorMadLoop):
                        relative_position_of_root_process_dir = '../..'
                        ME_name = '_%s'%me_accessor.proc_dir
                        source_dir = '$(PROC_ROOT)/%s/SubProcesses/' % me_accessor.proc_dir
                    else:
                        relative_position_of_root_process_dir = '../../..'
                        ME_name = '_%s__%s'%(me_accessor.proc_dir, me_accessor.proc_name)
                        source_dir = '$(PROC_ROOT)/%s/SubProcesses/P%s/' % (me_accessor.proc_dir, me_accessor.proc_name)
                    self.build_info['matrix_elements'].append({
                        'lib_name': me_accessor.get_library_name(),
                        'source_dir': source_dir,
                        'relative_position_of_root_process_dir' : relative_position_of_root_process_dir,
                        'makefile_env_variables': {
                            'MEDIR': '%s/lib/matrix_elements/'%relative_position_of_root_process_dir,
                            'MENAME': ME_name
                        },
                        'makefile_target': '%s/lib/matrix_elements/lib%s.a' % ( relative_position_of_root_process_dir,
                                                                                        me_accessor.get_library_name() )
                    })
                rust_writer.writelines(ME_rs_template, context={}, replace_dictionary=repl_dict)

            # Logic for rust export of low level subtraction currents
            if isinstance(me_accessor, accessors.SubtractionCurrentAccessor):
                # No pre-processing per subtraction counterterm to be performed at this stage
                pass

        rust_writer.close()

        return

    def translate_low_level_code(self, low_level_code, function_prefix=''):
        """ Translates the metacode for calling matrix element or computing counterterm into rust code."""

        rust_code = []
        necessary_imports = []
        replacement_dict = {'function_prefix':function_prefix}

        def translate_argument(arg):
            if isinstance(arg, tuple) and arg[0]=='EpsilonExpansion':
                return "EpsilonExpansion::from_slice(&[%s])"%(
                    ','.join('(%d, %s)'%(k, translate_argument(v)) for k, v in sorted(arg[1].items(), key=lambda k: k[0]))
                )
            elif isinstance(arg, list):
                return '[%s]'%(','.join(translate_argument(element) for element in arg))
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
                    necessary_imports.append(instruction[2])
            elif instruction[0] == 'set':
                # Assume immutable for now
                rust_code.append('let %s = %s;'%(instruction[1], translate_argument(instruction[2])))
            elif instruction[0] == 'call':
                rust_code.append('%s;'%translate_function_call(instruction[1],instruction[2]))
            elif instruction[0] == 'call_and_set':
                rust_code.append('let %s = %s;'%(instruction[1], translate_function_call(instruction[2],instruction[3])))
            elif instruction[0] == 'return':
                rust_code.append('%s'%translate_argument(instruction[1]))
            elif instruction[0] == 'print':
                rust_code.append('println!(%s);'%translate_argument(instruction[1]))

        return rust_code, necessary_imports

    def write_rust_counterterm(self, integrand, runtime_symbolic_inputs, i_process, process_key, i_CT, CT,
            counterterm_evaluators_rust_writer, current_evaluators_rust_writer, i_input_mapping=None, CT_input_mapping=None):
        """
        Returns
                headers, instantiation
        Where "headers" is a list of include statement necessary for that particular instantiation while
        "instantiation" is a block of code instantiating the specified local counterterm together with
        its evaluators.

        Given the two rust writer specified, this function writes the implementation of the evaluator
        for the counterterm specified as well as of all the evaluator of all the currents it contains.

        The options i_input_mapping and CT_input_mapping refer to the input mappings of the integrated counterterms
        and they are set to None for local counterterms.
        """

        counterterm_short_name = None
        if i_input_mapping is not None:
            counterterm_short_name += '%s_P%d_ICT%d_m%d' % (
                    self.get_integrand_short_name(integrand), i_process, i_CT,i_input_mapping)
        else:
            counterterm_short_name = '%s_P%d_CT%d' % (self.get_integrand_short_name(integrand), i_process, i_CT)

        counterterm_evaluator_context = {}
        counterterm_evaluator_repl_dict = {
            'counterterm_short_name' : counterterm_short_name,
            'integrand_id' : integrand.ID,
            'process_id' : i_process,
            'counterterm_id' : i_CT,
            'mapping_id' : 0 if i_input_mapping is None else i_input_mapping,
            'ME_imports' : None,
            'ME_calls' : None,
            'counterterm_evaluator_definition' : None,
            'counterterm_evaluator_construction': None,
        }

        counterterm_instantiation_includes = []
        counterterm_instantiation = [
            'Counterterm::new(%(integrand_id)d, %(process_id)d, %(CT_id)d, %(CT_mapping_id)d,'%{
                'integrand_id' : integrand.ID,
                'process_id' : i_process,
                'CT_id' : i_CT,
                'CT_mapping_id' : -1 if i_input_mapping is None else i_input_mapping
            },
            '&param_card_path','run_card', 'param_card', 'settings_card',
            'integrands:%sCountertermEvaluator(&param_card_path)'%counterterm_short_name,
        ]
        if i_input_mapping is None:
            header_lines_for_subtraction_currents, subtraction_currents_instantiation, \
            counterterm_evaluation_code, ME_calls_code =  integrand.evaluate_local_counterterm(
                CT, runtime_symbolic_inputs[i_process]['PS_point'],
                runtime_symbolic_inputs[i_process]['base_weight'],
                runtime_symbolic_inputs[i_process]['mu_r'],
                runtime_symbolic_inputs[i_process]['mu_f1'],
                runtime_symbolic_inputs[i_process]['mu_f2'],
                runtime_symbolic_inputs[i_process]['xb_1'],
                runtime_symbolic_inputs[i_process]['xb_2'],
                runtime_symbolic_inputs[i_process]['xi1'],
                runtime_symbolic_inputs[i_process]['xi2'],
                integrand.processes_map[process_key][0].get('has_mirror_process'),
                runtime_symbolic_inputs[i_process]['all_flavor_configurations'][0],
                hel_config=None, # MC over helicities not supported yet.
                # The following options will remain as arguments of the function
                sector='sector',
                apply_flavour_blind_cuts='apply_flavour_blind_cuts',
                boost_back_to_com='boost_back_to_com',
                always_generate_event='always_generate_event',
                low_level_generation=True,
                current_evaluators_writer=current_evaluators_rust_writer
            )
        else:
            # notice that in this case 'CT' is and 'integrated_CT_characteristics' dictionary
            header_lines_for_subtraction_currents, subtraction_currents_instantiation, \
            counterterm_evaluation_code, ME_calls_code = integrand.evaluate_integrated_counterterm(
                CT, runtime_symbolic_inputs[i_process]['PS_point'],
                runtime_symbolic_inputs[i_process]['base_weight'],
                runtime_symbolic_inputs[i_process]['mu_r'],
                runtime_symbolic_inputs[i_process]['mu_f1'],
                runtime_symbolic_inputs[i_process]['mu_f2'],
                runtime_symbolic_inputs[i_process]['xb_1'],
                runtime_symbolic_inputs[i_process]['xb_2'],
                runtime_symbolic_inputs[i_process]['xi1'],
                runtime_symbolic_inputs[i_process]['xi2'],
                runtime_symbolic_inputs[i_process]['xi2'],
                CT['input_mappings'][i_input_mapping],
                runtime_symbolic_inputs[i_process]['all_flavor_configurations'][0],
                # The following options will remain as arguments of the function
                sector='sector',
                hel_config=None,
                compute_poles='compute_poles',
                low_level_generation=True,
                current_evaluators_writer=current_evaluators_rust_writer
            )
        counterterm_evaluator_repl_dict['counterterm_evaluation_code'] = counterterm_evaluation_code

        counterterm_instantiation.append(subtraction_currents_instantiation)
        counterterm_instantiation_includes.extend(header_lines_for_subtraction_currents)


        # Now build the ME_calls placeholder from the instruction set returned by the evaluate_counterterm function
        ME_calls_low_level_code_for_CT_evaluator = []
        ME_imports_for_CT_evaluator = []
        ME_evaluators_definition = []
        ME_evaluators_construction = []
        for (i_call, matrix_element_name), call_instructions in ME_calls_code.items():
            ME_calls_low_level_code_for_CT_evaluator.append('(%d, %d, %d, %d) => {'%(
                i_process, i_CT, -1 if i_input_mapping is None else i_input_mapping, i_call)
            )
            rust_code, necessary_imports = self.translate_low_level_code(call_instructions)
            ME_calls_low_level_code_for_CT_evaluator.extend(rust_code)
            ME_imports_for_CT_evaluator.extend(necessary_imports)
            ME_evaluators_definition.append(
                '%(MEname)sMatrixElement_evaluator: MatrixElementEvaluator<%(MEname)sMatrixElement>'%{
                    'MEname' : matrix_element_name
                }
            )
            ME_evaluators_construction.append(
                '%(MEname)sMatrixElement_evaluator: MatrixElementEvaluator::new(card_filename, %(MEname)sMatrixElement {})'%{
                    'MEname': matrix_element_name
                }
            )

        counterterm_evaluator_repl_dict['ME_imports'] = '\n'.join(ME_imports_for_CT_evaluator)
        counterterm_evaluator_repl_dict['ME_calls'] = '\n'.join(ME_calls_low_level_code_for_CT_evaluator)
        counterterm_evaluator_repl_dict['counterterm_evaluator_definition'] = ',\n'.join(ME_evaluators_definition)
        counterterm_evaluator_repl_dict['counterterm_evaluator_construction'] = ',\n'.join(ME_evaluators_construction)


        counterterm_evaluators_rust_writer.writelines(
            open(pjoin(self.dynamic_template_path, 'counterterm_evaluator.rs'), 'r').read(),
            context=counterterm_evaluator_context, replace_dictionary=counterterm_evaluator_repl_dict
        )

        return counterterm_instantiation_includes, ','.join(counterterm_instantiation)

    def write_integrand_instantiator(self, integrand, runtime_symbolic_inputs, repl_dict):
        """ Export the integrand_instantiator.rs code for the particular integrand in argument."""

        integrand_short_name = self.get_integrand_short_name(integrand)
        integrand_export_path = pjoin(self.export_dir, 'rust', 'madnklo', 'src', 'integrands', integrand_short_name)

        integrand_instantiation = []
        integrand_instantiation_header = []
        integrand_instantiation.append("// Now instantiating integrand '%s'"%integrand_short_name)
        integrand_instantiation.append('integrands.insert(%d, Integrand::new('%integrand.ID)
        instantiation_repl_dict = {}
        instantiation_repl_dict['n_processes'] = len(integrand.processes_map)
        instantiation_repl_dict['n_initial'] = len(integrand.masses[0])
        instantiation_repl_dict['n_final'] = len(integrand.masses[1])
        representative_process = integrand.processes_map.values()[0][0]
        masses_symbols = [
            tuple('param_card.%s'%self.model.get_particle(pdg_code).get('mass') for pdg_code in representative_process.get_initial_ids()),
            tuple('param_card.%s'%self.model.get_particle(pdg_code).get('mass') for pdg_code in representative_process.get_final_ids()),
        ]
        instantiation_repl_dict['masses'] = "(vec![%s],vec![%s])"%(
            ','.join('%s'%m for m in masses_symbols[0]),','.join('%s'%m for m in masses_symbols[1]))
        instantiation_repl_dict['all_flavor_configurations'] = "hashmap![%s]"%(','.join('%d => %s'%
            (i_process, "vec![%s]"%(','.join('(vec![%s],vec![%s])'%(
                ','.join('%d'%pdg for pdg in flavor_config[0]),
                ','.join('%d'%pdg for pdg in flavor_config[1])) for flavor_config in symbolic_inputs['all_flavor_configurations'])))
            for i_process, symbolic_inputs in sorted(runtime_symbolic_inputs.items(), key = lambda el: el[0]) ) )
        instantiation_repl_dict['integrand_evaluator'] = "Box::new(IntegrandEvaluator_%s::new(&param_card_path))"%integrand_short_name

        all_local_CT_headers = []
        all_integrated_CT_headers = []
        instantiation_repl_dict['local_counterterms'] = 'vec!['
        instantiation_repl_dict['integrated_counterterms'] = 'vec!['

        # Now open the files that will be used for the evaluators of the counterterms and its subtraction currents.

        local_counterterm_evaluators_rust_writer = writers.RustWriter(
                                           pjoin(integrand_export_path, 'all_local_counterterm_evaluators.rs'), opt='w')
        local_counterterm_evaluators_rust_writer.writelines(
                                    open(pjoin(self.dynamic_template_path, 'all_counterterm_evaluators.rs'), 'r').read())
        integrated_counterterm_evaluators_rust_writer = writers.RustWriter(
                                           pjoin(integrand_export_path, 'all_integrated_counterterm_evaluators.rs'), opt='w')
        integrated_counterterm_evaluators_rust_writer.writelines(
                                    open(pjoin(self.dynamic_template_path, 'all_counterterm_evaluators.rs'), 'r').read())

        local_current_evaluators_rust_writer = writers.RustWriter(
                                            pjoin(integrand_export_path, 'all_local_subtraction_current_evaluators.rs'), opt='w')
        local_current_evaluators_rust_writer.write(
                        self.subtraction_module.exporter.rust_template_files['all_subtraction_current_evaluators.rs'])
        integrated_current_evaluators_rust_writer = writers.RustWriter(
                                            pjoin(integrand_export_path, 'all_integrated_subtraction_counterterm_evaluators.rs'), opt='w')
        integrated_current_evaluators_rust_writer.write(
                        self.subtraction_module.exporter.rust_template_files['all_subtraction_current_evaluators.rs'])


        # list of dictionary of process information that must be stored for each process contributing to that integrand
        process_info = []
        for i_process, (process_key, (process, mapped_processes)) in enumerate(sorted(integrand.processes_map.items())):

            # Save process information
            process_info.append({
                'n_unresolved_particles' : integrand.contribution_definition.n_unresolved_particles,
                'has_mirror_process' : process.get('has_mirror_process')
            })

            # Extract the call for the physical matrix element event

            # Now write local counterterms
            instantiation_repl_dict['local_counterterms'] += 'vec!['
            if integrand.has_local_counterterms():
                instantiation_repl_dict['local_counterterms'] += '\n\\ Local counterterms'
                local_counterterms_instantiation = []
                for i_CT, local_CT in enumerate(integrand.counterterms[process_key]):
                    # Write the source code for the evaluator of the local_CT and that of all its contributing subtraction current.
                    local_CT_headers, local_CT_instantiation = self.write_rust_counterterm(
                        integrand, runtime_symbolic_inputs, i_process, process_key, i_CT, local_CT,
                        local_counterterm_evaluators_rust_writer,
                        local_current_evaluators_rust_writer
                    )
                    all_local_CT_headers.extend(local_CT_headers)
                    local_counterterms_instantiation.append(local_CT_instantiation)

                instantiation_repl_dict['local_counterterms'] += ',\n'.join(local_counterterms_instantiation)
            else:
                instantiation_repl_dict['local_counterterms'] += 'vec![],'
            instantiation_repl_dict['local_counterterms'] += '],'

            # And integrated ones
            instantiation_repl_dict['integrated_counterterms'] = 'vec!['
            if integrand.has_integrated_counterterms():
                instantiation_repl_dict['integrated_counterterms'] += '\n\\ Integrated counterterms'
                instantiation_repl_dict['integrated_counterterms'] = 'vec!['
                for i_CT, counterterm_characteristics in enumerate(integrand.integrated_counterterms[process_key]):
                    integrated_counterterms_instantiation = []
                    for i_mapping, input_mapping in enumerate(counterterm_characteristics['input_mappings']):
                        integrated_CT_headers, integrated_CT_instantiation = self.write_rust_counterterm(
                            integrand, runtime_symbolic_inputs, i_process, process_key, i_CT, counterterm_characteristics,
                            integrated_counterterm_evaluators_rust_writer,
                            integrated_current_evaluators_rust_writer,
                            i_input_mapping = i_mapping, CT_input_mapping = input_mapping
                        )
                        integrated_counterterms_instantiation.append(integrated_CT_instantiation)
                        all_integrated_CT_headers.extend(integrated_CT_headers)
                    instantiation_repl_dict['integrated_counterterms'] += ',\n'.join(integrated_counterterms_instantiation)
                instantiation_repl_dict['integrated_counterterms'] += '],'
            else:
                instantiation_repl_dict['integrated_counterterms'] = 'vec![vec![]],'
            instantiation_repl_dict['integrated_counterterms'] += '],'

        # Close streams in which counterterms and subtraction current evaluators are written.
        local_counterterm_evaluators_rust_writer.close()
        integrated_counterterm_evaluators_rust_writer.close()
        local_current_evaluators_rust_writer.close()
        integrated_current_evaluators_rust_writer.close()

        instantiation_repl_dict['process_info'] = 'vec![' + ','.join('ProcessInfo::new(%d,%s)'%(
                        info['n_unresolved_particles'], 'true' if info['has_mirror_process'] else 'false'
                                                                                    ) for info in process_info ) + ']'
        integrand_instantiation.append(
            (   '%(n_processes)d, %(all_flavor_configurations)s, '+
                '%(process_info)s, %(n_initial)d, %(n_final)d, %(masses)s, '+
                'run_card, param_card, settings_card, %(integrand_evaluator)s,'+
                '\n%(local_counterterms)s,'+
                '\n%(integrated_counterterms)s'
            )%instantiation_repl_dict
        )

        integrand_instantiation.append(')')

        integrand_instantiation_header.append('// Integrand evaluator')
        integrand_instantiation_header.append('use crate::integrands::integrand_evaluator_%s::IntegrandEvaluator_%s;'%(integrand_short_name,integrand_short_name))
        integrand_instantiation_header.append('// Local counterterm evaluators')
        integrand_instantiation_header.extend(all_local_CT_headers)
        integrand_instantiation_header.append('// Integrated counterterm evaluator')
        integrand_instantiation_header.extend(all_integrated_CT_headers)

        # Now update the place holder
        repl_dict['instantiate_integrands'] += '\n'.join(integrand_instantiation) + ');' + '\n'
        repl_dict['instantiate_integrands_header'] += '\n'.join(integrand_instantiation_header) + '\n'

        repl_dict['integrands_include'] += 'pub mod integrand_evaluator_%s;\n' % integrand_short_name


    def write_integrand_evaluator(self, integrand, runtime_symbolic_inputs, repl_dict):
        """ Export the integrand evaluator code in in integrands/integrand_evaluator_B_1.rs"""

        integrand_short_name = self.get_integrand_short_name(integrand)
        integrand_export_path = pjoin(self.export_dir, 'rust', 'madnklo', 'src', 'integrands',integrand_short_name)

        integrand_evaluator_repl_dict = {}
        integrand_evaluator_repl_dict['integrand_ID'] = integrand.ID
        integrand_evaluator_repl_dict['integrand_short_name'] = integrand_short_name
        integrand_evaluator_repl_dict['integrand_evaluator_definition'] = ''
        integrand_evaluator_repl_dict['integrand_evaluator_construction'] = ''

        ME_calls_instructions = {}
        for i_process, (process_key, (process, mapped_processes)) in enumerate(sorted(integrand.processes_map.items())):

            call_key = (i_process,)

            # Extract the call for the physical matrix element event
            # Reconstruct what are the necessary ME calls in this sigma call so as to be able to hardcode them
            low_level_calling_code = integrand.all_MEAccessors(
                process,
                runtime_symbolic_inputs[i_process]['PS_point'],
                runtime_symbolic_inputs[i_process]['alpha_s'],
                runtime_symbolic_inputs[i_process]['mu_r'],
                pdgs=runtime_symbolic_inputs[i_process]['all_flavor_configurations'][0],
                low_level_code_generation=True)

            ME_calls_instructions[call_key] = low_level_calling_code

        ME_calls_lines = []
        evaluator_header = []
        evaluator_definition = []
        evaluator_construction = []
        matrix_element_imports = []
        for call_key, low_level_code in sorted(ME_calls_instructions.items(), key=lambda k: k[0]):
            ME_calls_lines.append('(%d,) => {'%call_key)
            rust_code, needed_matrix_elements = self.translate_low_level_code(low_level_code, function_prefix=repl_dict['C_binding_prefix'])
            for mat in needed_matrix_elements:
                matrix_element_imports.append('use crate::matrix_elements::%s;' % mat)
                evaluator_definition.append('%s_evaluator: MatrixElementEvaluator<%s>,' % (mat, mat))
                evaluator_construction.append('%s_evaluator: MatrixElementEvaluator::new(card_filename, %s {}),' % (mat, mat))
            ME_calls_lines.extend(rust_code)
            ME_calls_lines.append('},')
        integrand_evaluator_repl_dict['ME_calls'] = '\n'.join(ME_calls_lines)

        rust_writer = writers.RustWriter(pjoin(integrand_export_path,'integrand_evaluator_%s.rs'%integrand_short_name),opt='w')
        integrand_evaluator_context = {}
        integrand_evaluator_repl_dict['header'] = '\n'.join(evaluator_header)
        integrand_evaluator_repl_dict['ME_imports'] = '\n'.join(matrix_element_imports)
        integrand_evaluator_repl_dict['integrand_evaluator_definition'] = '\n'.join(evaluator_definition)
        integrand_evaluator_repl_dict['integrand_evaluator_construction'] = '\n'.join(evaluator_construction)
        rust_writer.writelines(
            open(pjoin(self.dynamic_template_path, 'integrand_evaluator.rs'), 'r').read(),
            context=integrand_evaluator_context, replace_dictionary=integrand_evaluator_repl_dict)
        rust_writer.close()

    def get_integrand_short_name(self, integrand):
        """ Helper function to centralise the choice of the integrand short name."""
        return '%s_%d' % (integrand.get_short_name(), integrand.ID)

    def export(self, integrand, repl_dict):
        """ Export one particular integrand. """
        #import pdb
        #pdb.set_trace()

        # Build the directory that will contain all rust resources pertaining to that integrand
        integrand_short_name = self.get_integrand_short_name(integrand)
        integrand_export_path = pjoin(self.export_dir, 'rust', 'madnklo', 'src', 'integrands',integrand_short_name)
        os.makedirs(integrand_export_path)

        # Make sure to first disable accessor caches
        accessors.deactivate_cache()

        # Build symbolic inputs to use when going through the low-level generation
        runtime_symbolic_inputs = {}
        runtime_symbolic_inputs_common_to_all_i_process = {
            'alpha_s'       : 'alpha_s',
            'mu_r'          : 'mu_r',
            'base_weight'   : 'base_weight',
            'mu_r'          : 'mu_r',
            'mu_f1'         : 'mu_f1',
            'mu_f2'         : 'mu_f2',
            'xb_1'          : 'xb_1',
            'xb_2'          : 'xb_2',
            'xi1'           : 'xi1',
            'xi2'           : 'xi2',
        }
        for i_process, (process_key, (process, mapped_processes)) in enumerate(sorted(integrand.processes_map.items())):

            # Build an empty dictionary for the symbolic inputs of the i_process
            runtime_symbolic_inputs[i_process] = dict(runtime_symbolic_inputs_common_to_all_i_process)

            # First construct the list of
            process_pdgs = process.get_cached_initial_final_pdgs()
            all_processes = [process,]+mapped_processes
            all_flavor_configurations = []
            # The process mirroring is accounted for at the very end only
            for proc in all_processes:
                initial_final_pdgs = proc.get_cached_initial_final_pdgs()
                all_flavor_configurations.append(initial_final_pdgs)
            runtime_symbolic_inputs[i_process]['all_flavor_configurations'] = all_flavor_configurations

            # We can now use the all_flavor_configurations built here to hard-code it in the rust template.
            PS_point = {}
            for i_leg, leg in enumerate(process.get_initial_legs()):
                PS_point[leg.get('number')] = 'p[&%d]'%leg.get('number')
            for i_leg, leg in enumerate(process.get_final_legs()):
                PS_point[leg.get('number')] = 'p[&%d]'%leg.get('number')
            runtime_symbolic_inputs[i_process]['PS_point'] = PS_point

        # We must now hard-code the collected low-level code into the dynamic rust templates
        # ==================================================================================

        # First output the code for instantiating this particular integrand in the AllIntegrands struct
        self.write_integrand_instantiator(integrand, runtime_symbolic_inputs, repl_dict)

        # Now we must write out the corresponding integrand evaluator
        self.write_integrand_evaluator(integrand, runtime_symbolic_inputs, repl_dict)

        return

    @staticmethod
    def compile(root_path):
        """ Compile the rust backend."""
        # Note: this should be done as much as possible using makefile targets so that one can easily recompile
        # the entire distribution manually.
        # So for now we can directly use the cmd 'make' from within the rust directory.
        return misc.compile(cwd=pjoin(root_path, 'rust'))

    @staticmethod
    def write_rust_settings_card(root_path, options):
        """ Write the rust settings card from MG5aMC options and the process root path."""

        try:
            import yaml
            from yaml import Loader, Dumper
            noalias_dumper = Dumper
            noalias_dumper.ignore_aliases = lambda self, data: True
        except ImportError:
            raise MadGraph5Error("The pyYAML python dependency is necessary for exporting %s in the '%s' format."%(
                                                                                       self.__class__.__name__, format))

        settings = {'root_path' : root_path}

        try:
            lhapdf_libdir = options['lhapdf'] # subprocess.Popen([options['lhapdf'], '--libdir'], stdout=subprocess.PIPE).stdout.read().strip()
        except:
            raise InvalidCmd("Could not determine the location of the LHAPDF library. Verify the value of the 'lhapdf' option of MG5aMC.")

        settings['lhapdf_library_path'] = os.path.abspath(pjoin(lhapdf_libdir,'libLHAPDF.a'))

        return yaml.dump(settings, Dumper=noalias_dumper, default_flow_style=False)

    def write_rust_makefile(self):
        """ Write the makefile for steering the compilation of all resources necessary to rust as well as the rust
        executable."""

        makefile_repl_dict = {}

        makefile_repl_dict['helas_libraries'] = '\\\n '.join(
            '$(PROC_ROOT)/lib/helas/lib%s.a'%info['lib_name'] for info in self.build_info['helas'])
        makefile_repl_dict['matrix_element_libraries'] = '\\\n '.join(
            '$(PROC_ROOT)/lib/matrix_elements/lib%s.a'%info['lib_name'] for info in self.build_info['matrix_elements'])


        helas_targets = []
        for info in self.build_info['helas']:
            helas_targets.append('$(PROC_ROOT)/lib/helas/lib%s.a: %s'%(info['lib_name'],info['source_dir']))
            helas_targets.append('\techo "Compiling %s ..."'%(info['source_dir'].replace('$(PROC_ROOT)','')))
            helas_targets.append('\tcd %s && %s make %s >>  ../../../rust/rust_compilation.log;'%(
                info['source_dir'], ' '.join('%s=%s'%(k,v) for k,v in info['makefile_env_variables'].items()),
                info['makefile_target']
            ))
        makefile_repl_dict['helas_targets'] = '\n'.join(helas_targets)

        matrix_element_targets = []
        for info in self.build_info['matrix_elements']:
            matrix_element_targets.append('$(PROC_ROOT)/lib/matrix_elements/lib%s.a: %s'%(info['lib_name'],info['source_dir']))
            matrix_element_targets.append('\techo "Compiling %s ..."'%(info['source_dir'].replace('$(PROC_ROOT)','')))
            matrix_element_targets.append('\tcd %s && %s make %s >>  %s/rust/rust_compilation.log;'%(
                info['source_dir'], ' '.join('%s=%s'%(k,v) for k,v in info['makefile_env_variables'].items()),
                info['makefile_target'], info['relative_position_of_root_process_dir']
            ))
        makefile_repl_dict['matrix_element_targets'] = '\n'.join(matrix_element_targets)

        clean_helas = []
        for info in self.build_info['helas']:
            clean_helas.append('\trm -f $(PROC_ROOT)/lib/helas/lib%s.a'%info['lib_name'])
            clean_helas.append('\tcd %s && make clean'%info['source_dir'])
        makefile_repl_dict['clean_helas'] = '\n'.join(clean_helas)

        clean_matrix_elements = []
        for info in self.build_info['matrix_elements']:
            clean_matrix_elements.append('\trm -f $(PROC_ROOT)/lib/matrix_elements/lib%s.a' % info['lib_name'])
            clean_matrix_elements.append('\tcd %s && make clean' % info['source_dir'])
        makefile_repl_dict['clean_matrix_elements'] = '\n'.join(clean_matrix_elements)

        template = open(pjoin(self.dynamic_template_path, 'Makefile.template'),'r').read()
        return template%makefile_repl_dict

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

        # write the module inclusion for all the integrands
        rust_writer = writers.RustWriter(pjoin(rust_export_path, 'integrands', 'mod.rs'), opt='w')
        rust_writer.writelines(
            r"%(integrands_include)s",
            context={}, replace_dictionary=repl_dict)
        rust_writer.close()

        # Also write the yaml equivalent of the cards, note that these will always be kept in sync before
        # launching in the ME7 command-line interface

        # The yaml-formatted param card
        open(pjoin(self.export_dir, 'rust','Cards','param_card.yaml'),'w').write(
            self.model.write_yaml_parameters_card())

        # The yaml-formatted run card
        open(pjoin(self.export_dir, 'rust','Cards','run_card.yaml'),'w').write(
            banner_mod.RunCardME7(pjoin(self.export_dir, 'Cards', 'run_card.dat')).export(
                pjoin(self.export_dir, 'Cards', 'run_card.dat'), format='yaml'))

        # The yaml-formatted settings card
        open(pjoin(self.export_dir, 'rust','Cards','settings.yaml'),'w').write(
            RustExporter.write_rust_settings_card(self.export_dir, self.options))

        # Now write the rust makefile for building the rust exectuable using cargo
        open(pjoin(self.export_dir, 'rust','Makefile'),'w').write(self.write_rust_makefile())