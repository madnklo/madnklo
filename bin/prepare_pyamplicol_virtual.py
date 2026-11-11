#!/usr/bin/env python3
"""Prepare an ordinary QCD virtual pilot, separately from tree artifacts.

This diagnostic command does not enable generated MadNkLO virtual dispatch.
Scheme and UV/R2 completeness are explicit caller declarations, not inferred
from a model name. Symbolic work stays inside pyAmpliCol.
"""
import argparse
import copy
from dataclasses import replace
from fractions import Fraction
import json
import itertools
from pathlib import Path
import subprocess

from prepare_pyamplicol import export_metadata, load_records, prepare_model, selected_sdk


def load_virtual_records(directory):
    records = {}
    for path in sorted(directory.glob('*/VirtualProvider/*.json')):
        record = json.loads(path.read_text())
        pdgs, ninitial = record.get('pdgs'), record.get('ninitial')
        parameters, requirements = record.get('external_parameters'), record.get('required_correlations')
        if (not isinstance(pdgs, list) or len(pdgs) < 3 or
                any(type(pdg) is not int or pdg == 0 for pdg in pdgs) or
                type(ninitial) is not int or not 0 < ninitial < len(pdgs) or
                record.get('loop_coupling_orders') != ['QCD'] or
                not isinstance(parameters, list) or not parameters or
                any(not isinstance(name, str) or not name for name in parameters) or
                len(set(parameters)) != len(parameters) or
                not isinstance(requirements, dict) or
                set(requirements) != {'color_order','spin_order'} or
                any(type(order) is not int or not 0 <= order <= 2 for order in requirements.values())):
            raise ValueError('Invalid virtual manifest: '+str(path))
        expected = 'p%d_%s' % (ninitial, '_'.join(map(str,pdgs)))
        if record.get('key') != expected:
            raise ValueError('Invalid virtual manifest: '+str(path))
        previous = records.get(expected)
        if previous is not None:
            if ({k:v for k,v in previous.items() if k != 'required_correlations'} !=
                    {k:v for k,v in record.items() if k != 'required_correlations'}):
                raise ValueError('Conflicting virtual definitions for '+expected)
            record['required_correlations'] = {name:max(value,previous['required_correlations'][name])
                                              for name,value in record['required_correlations'].items()}
        records[expected] = record
    return records


def virtual_header(records, inspections, parameters, path):
    """Compile the exported Fortran parameter order into SDK name bindings."""
    lines, definitions = [], []
    for index, (key, record) in enumerate(sorted(records.items())):
        source = record['external_parameters']
        active = set(inspections[key]['parameter_names'])
        if active-set(source):
            raise ValueError('Virtual SDK parameters absent from generated Fortran: '+', '.join(sorted(active-set(source))))
        lines.append('static const VirtualParameter virtual_parameters_%d[] = {' % index)
        for position, name in enumerate(source):
            real, imaginary = parameters[name]
            if float(imaginary) != 0:
                raise ValueError('Generated virtual routing requires real parameters: '+name)
            lines.append('{%s,%d,%s,%.17g},' %
                         (json.dumps(name), position, 'true' if name in active else 'false', float(real)))
        lines.append('};')
        aliases = []
        catalogue = inspections[key].get('color_correlations',[])
        ids = [entry['id'] for entry in catalogue]
        proof = path.parent/(key+'.conventions.json')
        if record.get('legacy_colors'):
            if not proof.is_file():
                raise ValueError('Missing Born-loop colour compatibility proof for '+key)
            for row, alias in enumerate(json.loads(proof.read_text())['aliases']):
                terms = alias['terms']
                name = 'virtual_terms_%d_%d'%(index,row)
                lines.append('static const VirtualTerm %s[] = {'%name)
                lines += ['{%d,%.17g,%.17g},'%(ids.index(term['id']),float(Fraction(term['real'])),
                                                  float(Fraction(term['imag']))) for term in terms]
                if not terms:
                    lines.append('{0,0,0},')
                lines.append('};')
                aliases.append('{%d,%d,%d,%s}'%(alias['first'],alias['second'],len(terms),name))
        name = 'virtual_aliases_%d'%index
        lines += ['static const VirtualAlias %s[] = {'%name,',\n'.join(aliases) if aliases else '{0,0,0,nullptr}','};']
        definitions.append('{%s,%s,%d,virtual_parameters_%d,%d,%s,%d,%d,%s,%s}' %
            (json.dumps(key),json.dumps(inspections[key]['process_id']),len(source),index,
             len(aliases),name,len(ids),ids.index('born') if ids else 0,
             json.dumps(json.dumps(ids,separators=(',',':'))) if ids else 'nullptr',
             json.dumps(json.dumps(inspections[key].get('spin_correlations',[]),separators=(',',':'))) if ids else 'nullptr'))
    lines += ['const VirtualDefinition *virtual_definition(const std::string &key) {',
              'static const VirtualDefinition definitions[] = {', ',\n'.join(definitions), '};',
              'for (const auto &definition : definitions) if (key == definition.key) return &definition;',
              'return nullptr;', '}']
    path.write_text('\n'.join(lines)+'\n')


def loop_configuration(record, declaration):
    from pyamplicol.config import resolve_config

    config = resolve_config(dict(
        action='generate',
        process=dict(loop_coupling_orders=['QCD'], coupling_order_policy='explicit',
                     max_coupling_orders={name: value for name, value in record.get('born_orders',record['orders']).items()
                                          if name != 'WEIGHTED'}),
        loop=dict(r2_strategy='model', model=declaration),
        color=dict(accuracy='full', contraction='direct'),
        evaluator=dict(execution_mode='auto', backend='jit'),
        generation=dict(emit_api_bundle=False)))
    model = config.effective.loop.model.to_declaration()
    if (model is None or model.loop_coupling_orders != ('QCD',) or
            not model.uv_complete or not model.model_r2_complete):
        raise ValueError('The virtual pilot requires an explicit complete QCD UV/R2 declaration')
    return config


def virtual_parameters(parameters, inspection):
    """Use the loop source's parameter catalogue, not the empty tree-style view."""
    names = inspection['parameter_names']
    missing = set(names)-set(parameters)
    if missing:
        raise ValueError('Virtual model parameters missing from card: '+', '.join(sorted(missing)))
    for name, fixed in dict(inspection['fixed_parameters']).items():
        if name not in parameters or name in names:
            continue
        pair = fixed if isinstance(fixed, (tuple, list)) else [fixed, '0']
        if [Fraction(str(value)) for value in parameters[name]] != [Fraction(str(value)) for value in pair]:
            raise ValueError('Virtual fixed parameter differs from card: '+name)
    return {name: [str(value) for value in parameters[name]] for name in names}


def validate_virtual_orders(record, inspection):
    if 'virtual_squared_order' not in record:
        return
    expected = {name:value for name,value in record['virtual_squared_order'].items() if value}
    born = {name:2*value for name,value in record['born_orders'].items() if value}
    actual = inspection['orders']
    if (len(actual['born']) != 1 or len(actual['virtual']) != 1 or
            dict(actual['virtual'][0]['order']) != expected or dict(actual['born'][0]['order']) != born):
        raise ValueError('Prepared virtual coupling orders differ from the exported request')


def prepare(args, record_override=None, metadata_override=None, build_library=True):
    from pyamplicol import CorrelatorConfig, Generator, ModelSource, Runtime
    from pyamplicol.models.loading import compile_model_source
    from rusticol_config import load_sdk_info

    directory = args.directory.resolve()
    records = load_records(directory)
    if args.key not in records:
        raise ValueError('Unknown exported process key: '+args.key)
    record = record_override if record_override is not None else load_virtual_records(directory).get(args.key,records[args.key])
    declaration = json.loads(args.declaration.read_text())
    config = loop_configuration(record, declaration)
    sdk = selected_sdk(load_sdk_info, 'shared')
    metadata = metadata_override or (json.loads(args.metadata.read_text()) if getattr(args,'metadata',None)
                                     else export_metadata(directory, args.madnklo_python))
    root = directory/'VirtualProvider'
    root.mkdir(exist_ok=True)
    correlations = None
    if 'legacy_colors' in record and any(record['required_correlations'].values()):
        manifest = root/(args.key+'.conventions-input.json')
        proof = root/(args.key+'.conventions.json')
        manifest.write_text(json.dumps(record))
        subprocess.run([args.madnklo_python,str(Path(__file__).with_name('pyamplicol_virtual_legacy.py')),
                        str(manifest),str(proof)],check=True)
        conventions = json.loads(proof.read_text())
        gluons = [i+1 for i,pdg in enumerate(record['pdgs']) if pdg == 21]
        spins = [group for size in range(1,min(record['required_correlations']['spin_order'],len(gluons))+1)
                 for group in itertools.combinations(gluons,size)]
        correlations = dict(color_correlations=[entry for entry in conventions['catalogue'] if entry['bra']],
                            spin_correlations=spins)
    parameter_file = root/'parameters.json'
    if parameter_file.exists() and json.loads(parameter_file.read_text()) != metadata['parameters']:
        raise ValueError('Virtual parameters changed: use a fresh export')
    request = dict(record=record, loop_model=declaration, correlations=correlations,
                   loop_sm_scalar_compatibility=args.loop_sm_scalar_compatibility)
    saved = root/(args.key+'.record.json')
    artifact = root/args.key
    if artifact.exists():
        if not args.resume or not saved.is_file() or json.loads(saved.read_text()) != request:
            raise ValueError('Virtual artifact is not a resumable match: '+str(artifact))
    else:
        source = prepare_model(Path(metadata['model']), root, args.madnklo_python, include_loops=True)
        (source/'restrict_madnklo.dat').write_text(metadata['restriction_card'])
        compiled = compile_model_source(source, restriction='madnklo',
            cache_dir=root/'model-cache', require_supported=False)
        if args.loop_sm_scalar_compatibility:
            from pyamplicol.loops.model_compatibility import prepare_loop_sm_scalar_declarations
            prepared = prepare_loop_sm_scalar_declarations(compiled.raw_loop_model)
            compiled = replace(compiled, raw_loop_model=prepared.derived,
                               _serialized_path=None, _prepared_bundle=None)
        model_path = compiled.write(root/'loop-model.pyAmplicol-model.json')
        print('Generating ordinary virtual:', record['process'], flush=True)
        Generator(config).generate(record['process'], artifact,
            model=ModelSource.from_path(model_path),
            correlators=CorrelatorConfig.from_json_dict(correlations) if correlations else None)
    runtime = Runtime.load(artifact)
    inspection = runtime.inspect()
    validate_virtual_orders(record,inspection)
    exact_parameters = virtual_parameters(metadata['parameters'], runtime.inspect())
    runtime.set_model_parameters(exact_parameters)
    (root/(args.key+'.parameters.json')).write_text(json.dumps(exact_parameters, indent=2)+'\n')
    (root/(args.key+'.inspection.json')).write_text(json.dumps(dict(runtime.inspect()), indent=2)+'\n')
    saved.write_text(json.dumps(request, indent=2, sort_keys=True)+'\n')
    parameter_file.write_text(json.dumps(metadata['parameters'], indent=2, sort_keys=True)+'\n')
    if not build_library:
        return dict(runtime.inspect())
    library = root/'libmadnklo_pyamp_virtual.so'
    source_root = Path(__file__).resolve().parents[1]
    subprocess.run([args.cxx, '-std=c++17', '-O2', '-fPIC', '-shared',
        str(source_root/'madgraph/iolibs/pyamplicol_virtual.cc'),
        *sdk.cflags, *sdk.link_flags, '-o', str(library)], check=True)
    print('Virtual pilot artifact:', artifact)
    print('Virtual pilot library:', library)
    print('Generated MadNkLO V/RV routing is not enabled by this command.')


def prepare_virtual_artifacts(args, records, metadata):
    inspections = {}
    for key, record in sorted(records.items()):
        options = copy.copy(args)
        options.key = key
        options.declaration = args.loop_declaration
        inspections[key] = prepare(options, record, metadata, build_library=False)
    header = args.directory.resolve()/'VirtualProvider'/'virtual_parameters.inc'
    virtual_header(records, inspections, metadata['parameters'], header)
    return header


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    parser.add_argument('--key', required=True, help='exported BornProvider process key')
    parser.add_argument('--madnklo-python', required=True)
    parser.add_argument('--declaration', required=True, type=Path,
                        help='explicit loop-model scheme and UV/R2 completeness JSON')
    parser.add_argument('--loop-sm-scalar-compatibility', action='store_true',
                        help='apply pyAmpliCol\'s checked loop_sm scalar-arity preparation')
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--metadata', type=Path,
                        help='Diagnostic amplitude-only export metadata, not a full MadEvent export')
    parser.add_argument('--cxx', default='g++')
    args = parser.parse_args()
    try:
        prepare(args)
    except (ImportError, ValueError, RuntimeError, OSError, subprocess.CalledProcessError) as error:
        parser.exit(1, 'pyAmpliCol virtual preparation failed: %s\n'%error)


if __name__ == '__main__':
    main()
