#!/usr/bin/env python3
"""Prepare native tree artifacts for a Torino export (Python 3.11+)."""
import argparse
import ast
import json
import os
from pathlib import Path
import subprocess
import sys
import shutil
import itertools
import hashlib
import tempfile
import copy
from fractions import Fraction

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pyamplicol_conventions import dipole_product, outgoing_representations, validate_correlator


def prepare_model(source, root, python2, *, include_loops=False):
    destination = root / ('ufo_model_loop' if include_loops else 'ufo_model')
    if not destination.exists():
        shutil.copytree(source, destination,
                        ignore=shutil.ignore_patterns('*.pyc', '*.pkl', '__pycache__'))
        subprocess.run([python2, '-m', 'lib2to3', '-w', '-n', str(destination)],
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    if include_loops:
        return destination
    initializer = destination / '__init__.py'
    tree = ast.parse(initializer.read_text())
    class TreeOnly(ast.NodeTransformer):
        def visit_ImportFrom(self, node):
            node.names = [name for name in node.names if not name.name.startswith('CT_')]
            return node if node.names else None

        def visit_Assign(self, node):
            if any(isinstance(target, ast.Name) and target.id.startswith('all_CT')
                   for target in node.targets):
                return None
            return node
    tree = TreeOnly().visit(tree)
    initializer.write_text(ast.unparse(tree) + '\n')
    return destination


def load_records(directory):
    records = {}
    for path in sorted(directory.glob('*/BornProvider/*.json')):
        record = json.loads(path.read_text())
        record['orders'] = {name: value for name, value in record['orders'].items()
                            if name != 'WEIGHTED'}
        requirements = record.get('required_correlations')
        if (not isinstance(requirements, dict) or
                set(requirements) != {'color_order', 'spin_order'} or
                any(type(value) is not int or value not in (0, 1, 2)
                    for value in requirements.values())):
            raise ValueError('Missing or invalid subtraction requirements: regenerate '+str(path))
        previous = records.get(record['key'])
        if previous is None:
            records[record['key']] = record
            continue
        if ({key: value for key, value in previous.items() if key != 'required_correlations'} !=
                {key: value for key, value in record.items() if key != 'required_correlations'}):
            raise ValueError('Conflicting definitions for ' + record['key'])
        previous['required_correlations'] = {
            name: max(previous['required_correlations'][name], value)
            for name, value in requirements.items()}
    if not records:
        raise ValueError('No provider manifest: regenerate with MADNKLO_BORN_PROVIDER=pyamplicol')
    return records


def correlation_orders(record, args):
    """Manifest requirements are minima; explicit diagnostic flags may extend them."""
    required = record['required_correlations']
    legacy = max(required['color_order'], args.legacy_order or 0)
    return (max(required['color_order'], args.color_order or 0, legacy),
            max(required['spin_order'], args.spin_order or 0), legacy)


def load_color_requests(path, records):
    if path is None:
        return {}
    requests = json.loads(path.read_text())
    if not isinstance(requests, dict) or set(requests) - set(records):
        raise ValueError('Colour requests must map known process keys to declarations')
    for value in requests.values():
        if not isinstance(value, dict) or set(value)-{'color_correlations','spin_correlations'}:
            raise ValueError('Additional declarations may contain colour and spin correlations only')
    return requests


def color_legs(record):
    """Use model representations, not a PDG-number guess about coloured legs."""
    return sorted(outgoing_representations(record))


def quadruple_color_requests(record):
    """Both orders of every pair of full colour dipoles, including shared legs."""
    pairs = list(itertools.combinations_with_replacement(color_legs(record), 2))
    return [dipole_product(record, (first, second))
            for first in pairs for second in pairs]


def planned_catalogue(record, declarations):
    """Resolve public correlator declarations before compiling any amplitudes."""
    from pyamplicol.color.connections import ColorLeg
    from pyamplicol.color.connection_catalogue import build_color_correlator_catalogue

    operators = declarations.color_correlations
    if declarations.all_color_through_order is not None:
        outgoing = outgoing_representations(record)
        legs = tuple(ColorLeg(i+1, outgoing.get(i+1, 1)) for i in range(len(record['pdgs'])))
        operators += build_color_correlator_catalogue(legs,
            through_order=declarations.all_color_through_order)
    catalogue = [dict(id='born', bra=[], ket=[])] + [entry.to_json_dict() for entry in operators]
    for request in catalogue:
        validate_correlator(record, request)
    return catalogue


def catalogue_by_id(catalogue):
    result = {}
    for request in catalogue:
        identifier = request['id']
        if identifier in result:
            raise ValueError('Repeated colour correlation ID: '+identifier)
        result[identifier] = request
    return result


def legacy_catalogue(record, catalogue, proof):
    additional = proof.get('additional_requests', [])
    if not isinstance(additional, list):
        raise ValueError('Invalid additional legacy colour declarations')
    result = catalogue+additional
    for request in additional:
        validate_correlator(record, request)
    catalogue_by_id(result)
    return result


def legacy_header(maps,path):
    terms, aliases = [], []
    for key,value in sorted(maps.items()):
        for identifier,alias in sorted(value['aliases'].items()):
            aliases.append('{%s,%s,%d,%d}' % (json.dumps(key),json.dumps(identifier),len(terms),len(alias['terms'])))
            terms.extend('{%s,%.17g,%.17g}' % (json.dumps(term['id']),float(Fraction(term['real'])),
                                            float(Fraction(term['imag']))) for term in alias['terms'])
    # Constant arrays avoid thousands of optimized map/vector constructors.
    lines = ['const std::vector<LegacyTerm> *legacy_terms(const std::string &process,const std::string &id) {',
             'static const LegacyTerm terms[] = {', ',\n'.join(terms or ['{"",0,0}']), '};',
             'struct Alias { const char *process; const char *id; size_t offset; size_t count; };',
             'static const Alias aliases[] = {',
             ',\n'.join(aliases or ['{"","",0,0}']), '};',
             'size_t first=0,last=%d;' % len(aliases),
             'while(first<last) {',
             'size_t middle=first+(last-first)/2;',
             'int order=process.compare(aliases[middle].process);',
             'if(order==0) order=id.compare(aliases[middle].id);',
             'if(order>0) first=middle+1; else last=middle;', '}',
             'if(first==%d || process!=aliases[first].process || id!=aliases[first].id) return nullptr;' % len(aliases),
             'static std::map<size_t,std::vector<LegacyTerm>> cache;',
             'auto found=cache.find(first);',
             'if(found==cache.end()) {',
             'const auto &alias=aliases[first];',
             'found=cache.emplace(first,std::vector<LegacyTerm>(terms+alias.offset,terms+alias.offset+alias.count)).first;',
             '}', 'return &found->second;', '}']
    path.write_text('\n'.join(lines)+'\n')


def project_parameters(parameters, declarations):
    projected = {}
    for entry in declarations:
        if entry.name not in parameters:
            continue
        value = parameters[entry.name]
        if entry.mutable:
            projected[entry.name] = value
        elif complex(*value) != complex(entry.default_real,entry.default_imaginary):
            raise ValueError('Artifact fixed parameter differs from card: '+entry.name)
    return projected


def legacy_fingerprint(payload, source_root):
    """Invalidate exact proofs when their inputs or Python algebra change."""
    sources = [source_root/'bin'/name for name in
               ('prepare_pyamplicol.py', 'pyamplicol_legacy.py', 'pyamplicol_conventions.py')]
    sources += sorted((source_root/'madgraph').rglob('*.py'))
    digest = hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(',', ':')).encode())
    for path in sources:
        digest.update(str(path.relative_to(source_root)).encode()+b'\0')
        digest.update(hashlib.sha256(path.read_bytes()).digest())
    return digest.hexdigest()


def prepare_legacy(record, catalogue, args, root, source_root):
    payload = dict(record=record, catalogue=catalogue, order=args.legacy_order,
                   extension_limit=getattr(args, 'legacy_extension_limit', 1024))
    payload_path = root/(record['key']+'.legacy-input.json')
    destination = root/(record['key']+'.legacy.json')
    cache_path = root/(record['key']+'.legacy-cache.json')
    fingerprint = legacy_fingerprint(payload, source_root)
    if args.resume and destination.is_file() and cache_path.is_file():
        try:
            cache = json.loads(cache_path.read_text())
            serialized = destination.read_bytes()
            result = json.loads(serialized)
            if (cache.get('fingerprint') == fingerprint and
                    cache.get('result_sha256') == hashlib.sha256(serialized).hexdigest() and
                    result.get('process') == record['key']):
                payload_path.write_text(json.dumps(dict(payload,
                    catalogue=legacy_catalogue(record, catalogue, result))))
                print('Reusing exact legacy proof', record['key'], flush=True)
                return result
        except (OSError, ValueError, AttributeError):
            pass
    # Publish only a complete successful proof; a failed attempt keeps prior evidence.
    with tempfile.TemporaryDirectory(prefix='.legacy-', dir=root) as scratch:
        temporary = Path(scratch)/'proof.json'
        temporary_input = Path(scratch)/'input.json'
        temporary_input.write_text(json.dumps(payload))
        progress = getattr(args, 'legacy_progress', False)
        command = subprocess.run([args.madnklo_python, str(source_root/'bin/pyamplicol_legacy.py'),
                                  str(temporary_input), str(temporary)],
                                 stdout=None if progress else subprocess.PIPE,
                                 stderr=subprocess.STDOUT, text=True,
                                 env=dict(os.environ, MADNKLO_LEGACY_PROGRESS='1' if progress else '0'))
        if command.returncode:
            raise ValueError('Legacy colour matching failed:\n'+(command.stdout or 'See the streamed proof log above.'))
        serialized = temporary.read_bytes()
        result = json.loads(serialized)
        temporary_input.write_text(json.dumps(dict(payload,
            catalogue=legacy_catalogue(record, catalogue, result))))
        temporary_input.replace(payload_path)
        temporary.replace(destination)
        stamp = Path(scratch)/'cache.json'
        stamp.write_text(json.dumps(dict(fingerprint=fingerprint,
            result_sha256=hashlib.sha256(serialized).hexdigest()), indent=2)+'\n')
        stamp.replace(cache_path)
    return result


def selected_sdk(load_sdk_info, linkage):
    if linkage == 'static':
        return load_sdk_info()
    if linkage != 'shared':
        raise ValueError('Unknown SDK linkage: '+str(linkage))
    try:
        sdk = load_sdk_info(linkage='shared')
    except TypeError as error:
        raise ValueError('Install a pyAmpliCol SDK with shared-library support') from error
    if sdk.library.suffix not in ('.so', '.dylib'):
        raise ValueError('Shared SDK selection did not return a shared library')
    return sdk


def export_metadata(directory, python2):
    source_root = Path(__file__).resolve().parents[1]
    output = subprocess.check_output([
        python2, str(source_root/'bin/export_pyamplicol_parameters.py'),
        str(directory)], text=True)
    return json.loads(next(line.split(':', 1)[1] for line in output.splitlines()
                           if line.startswith('PYAMP_METADATA:')))


def prepare(args):
    from pyamplicol import ColorCorrelator, CorrelatorConfig, Generator, ModelSource, Runtime
    from pyamplicol.config import ColorConfig, EvaluatorConfig, ProcessConfig, RunConfig
    from rusticol_config import load_sdk_info

    if getattr(args,'native_rescue',False) and not hasattr(Runtime,'compile_correlated_rescue'):
        raise ValueError('Install the updated pyAmpliCol SDK before requesting native rescue')
    directory = args.directory.resolve()
    from prepare_pyamplicol_virtual import load_virtual_records, prepare_virtual_artifacts
    virtual_records = load_virtual_records(directory)
    linkage = getattr(args,'sdk_linkage',None) or ('shared' if virtual_records else 'static')
    if virtual_records and (linkage != 'shared' or not getattr(args,'loop_declaration',None)):
        raise ValueError('Combined preparation requires shared SDK linkage and --loop-declaration')
    sdk = selected_sdk(load_sdk_info, linkage)
    records = load_records(directory)
    requests = load_color_requests(args.color_requests, records)
    # Fail before generating artifacts or changing parameter snapshots.
    for key, record in records.items():
        color_legs(record)
        for request in requests.get(key, {}).get('color_correlations', []):
            validate_correlator(record, request)
    extra = {key: CorrelatorConfig.from_json_dict(value) for key, value in requests.items()}
    source_root = Path(__file__).resolve().parents[1]
    metadata = export_metadata(directory, args.madnklo_python)
    root = directory / 'BornProvider'
    root.mkdir(exist_ok=True)
    parameter_file = root / 'parameters.json'
    parameters = metadata['parameters']
    if parameter_file.exists() and json.loads(parameter_file.read_text()) != parameters:
        raise ValueError('Parameters changed: use a fresh export to avoid stale artifacts')
    parameter_file.write_text(json.dumps(parameters, indent=2, sort_keys=True) + '\n')
    model_path = prepare_model(Path(metadata['model']), root, args.madnklo_python)
    (model_path/'restrict_madnklo.dat').write_text(metadata['restriction_card'])
    model = ModelSource(kind='ufo', path=model_path, restriction='madnklo')
    legacy_maps = {}
    for record in records.values():
        color_order, spin_order, legacy_order = correlation_orders(record, args)
        artifact = root / record['key']
        colored = color_legs(record)
        gluons = [i + 1 for i, pdg in enumerate(record['pdgs']) if pdg == 21]
        operators = tuple(ColorCorrelator.dipole('T_%d_%d' % (a, b), a, b)
                          for a in colored for b in colored) if color_order else ()
        if color_order >= 2:
            operators += CorrelatorConfig.from_json_dict(dict(
                color_correlations=quadruple_color_requests(record))).color_correlations
        if record['key'] in extra:
            operators += extra[record['key']].color_correlations
        spins = tuple(group for size in range(1,min(spin_order,len(gluons))+1)
                      for group in itertools.combinations(gluons,size))
        if record['key'] in extra:
            spins = tuple(sorted(set(spins+extra[record['key']].spin_correlations)))
        declarations = CorrelatorConfig(color_correlations=operators,
            spin_correlations=spins,
            all_color_through_order=color_order if color_order>1 else None)
        correlated = bool(operators or spins or declarations.all_color_through_order)
        expected_catalogue = planned_catalogue(record, declarations) if correlated else []
        if legacy_order:
            print('Checking exact legacy compatibility for', record['process'], flush=True)
            legacy_args = copy.copy(args)
            legacy_args.legacy_order = legacy_order
            legacy_maps[record['key']] = prepare_legacy(record, expected_catalogue, legacy_args, root, source_root)
            proof = legacy_maps[record['key']]
            additions = proof.get('additional_requests', [])
            if additions:
                added = CorrelatorConfig.from_json_dict(dict(color_correlations=additions))
                declarations = CorrelatorConfig(
                    color_correlations=declarations.color_correlations+added.color_correlations,
                    spin_correlations=declarations.spin_correlations,
                    all_color_through_order=declarations.all_color_through_order)
                expected_catalogue = planned_catalogue(record, declarations)
                print('Declaring', len(additions), 'extra physical operators for exact legacy compatibility', flush=True)
            print('Matched', len(proof['aliases']), 'legacy IDs exactly', flush=True)
        orders = {name: value for name, value in record['orders'].items()
                  if name not in ('WEIGHTED',)}
        config = RunConfig(action='generate', color=ColorConfig(accuracy='full', contraction='direct'),
                           evaluator=EvaluatorConfig(execution_mode='compiled'),
                           process=ProcessConfig(max_coupling_orders=orders))
        print('Generating', record['process'], flush=True)
        print('  correlations: colour %d, spin %d, legacy %d' %
              (color_order, spin_order, legacy_order), flush=True)
        declaration_json = declarations.to_json_dict() if correlated else None
        sidecar = root/(record['key']+'.correlators.json')
        record_file = root/(record['key']+'.record.json')
        if artifact.exists():
            previous_record = json.loads(record_file.read_text()) if record_file.exists() else None
            legacy_input = root/(record['key']+'.legacy-input.json')
            if previous_record is None and legacy_input.exists():
                previous_record = json.loads(legacy_input.read_text())['record']
            if not args.resume or previous_record != record or not sidecar.exists() or \
                    json.loads(sidecar.read_text()) != declaration_json:
                raise ValueError('Existing artifact is not a verified resumable match: '+str(artifact))
            Runtime.load(artifact)
            print('Reusing verified artifact',record['key'],flush=True)
        else:
            Generator(config).generate(record['process'], artifact, model=model,
                correlators=declarations if correlated else None)
        record_file.write_text(json.dumps(record,indent=2,sort_keys=True)+'\n')
        (root / (record['key']+'.correlators.json')).write_text(
            json.dumps(declaration_json, indent=2, sort_keys=True)+'\n')
        runtime = Runtime.load(artifact)
        projected_parameters = project_parameters(parameters,runtime.physics.model_parameters)
        runtime.set_model_parameters({name:complex(*value) for name,value in projected_parameters.items()})
        (root/(record['key']+'.parameters.json')).write_text(
            json.dumps(projected_parameters,indent=2,sort_keys=True)+'\n')
        catalogue = [entry.to_json_dict() for entry in runtime.available_color_correlations()] if correlated else []
        if catalogue_by_id(catalogue) != catalogue_by_id(expected_catalogue):
            raise ValueError('Runtime colour catalogue differs from the exact preparation plan')
        (root/(record['key']+'.catalogue.json')).write_text(json.dumps(catalogue,indent=2,sort_keys=True)+'\n')
        if getattr(args,'native_rescue',False) and correlated:
            runtime.compile_correlated_rescue(root/(record['key']+'.rescue.so'),compiler=args.cxx)
    library = root / 'libmadnklo_pyamp.so'
    extra_flags = []
    if getattr(args,'native_rescue',False):
        extra_flags.append('-DMADNKLO_NATIVE_RESCUE')
    if legacy_maps:
        header = root/'legacy_correlations.inc'
        legacy_header(legacy_maps,header)
        extra_flags.append('-DMADNKLO_LEGACY_HEADER='+json.dumps(str(header)))
    virtual_sources = []
    if virtual_records:
        header = prepare_virtual_artifacts(args, virtual_records, metadata)
        extra_flags.append('-DMADNKLO_VIRTUAL_HEADER='+json.dumps(str(header)))
        virtual_sources = [str(source_root/'madgraph/iolibs'/name) for name in
                           ('pyamplicol_virtual.cc','pyamplicol_virtual_dispatch.cc')]
    subprocess.run([args.cxx, '-std=c++17', '-O2', '-fPIC', '-shared',*extra_flags,
        str(source_root / 'madgraph/iolibs/pyamplicol_native.cc'),
        *virtual_sources,
        *sdk.cflags, *sdk.link_flags, '-o', str(library)], check=True)
    print('Prepared provider. Set these variables before running the native executables:')
    import shlex
    print('export MADNKLO_PYAMP_ROOT=' + shlex.quote(str(root)))
    print('export MADNKLO_PYAMP_LIBRARY=' + shlex.quote(str(library)))
    if virtual_records:
        print('export MADNKLO_PYAMP_VIRTUAL_ROOT=' + shlex.quote(str(directory/'VirtualProvider')))
        print('Prepared ordinary and subtraction-required correlated virtuals.')
    if getattr(args,'native_rescue',False):
        print('export MADNKLO_PYAMP_CORRELATION_PRECISION=binary128')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    parser.add_argument('--madnklo-python', required=True, help='Absolute Python 2 interpreter path')
    parser.add_argument('--cxx', default='g++')
    parser.add_argument('--sdk-linkage', choices=('static','shared'),
                        help='SDK linkage (legacy tree preparation defaults to static)')
    parser.add_argument('--loop-declaration', type=Path,
                        help='Explicit UV/R2 and scheme declaration for exported virtuals')
    parser.add_argument('--loop-sm-scalar-compatibility', action='store_true',
                        help='Apply pyAmpliCol\'s checked loop_sm scalar-arity preparation')
    parser.add_argument('--native-rescue',action='store_true',
                        help='Compile optional Symbolica-free binary128 current/amplitude companions')
    parser.add_argument('--color-requests', type=Path,
                        help='Additional per-process pyAmpliCol colour declarations (JSON)')
    parser.add_argument('--color-order', type=int, choices=(0,1,2),
                        help='Diagnostic minimum colour order (default: contribution requirements)')
    parser.add_argument('--spin-order', type=int, choices=(0,1,2),
                        help='Diagnostic minimum spin order (default: contribution requirements)')
    parser.add_argument('--legacy-order',type=int,choices=(0,1,2),
                        help='Diagnostic minimum legacy-ID order (default: contribution requirements)')
    parser.add_argument('--legacy-extension-limit', type=int, default=1024,
                        help='Maximum extra ordered dipole products to test for exact legacy matching (0 disables)')
    parser.add_argument('--legacy-progress', action='store_true',
                        help='Stream exact colour-table and legacy-matching progress')
    parser.add_argument('--nnlo',action='store_true',
                        help='Compatibility flag; NNLO requirements are now read per contribution')
    parser.add_argument('--resume',action='store_true',
                        help='Reuse completed artifacts only when record, parameters and declarations match exactly')
    args = parser.parse_args()
    if args.legacy_extension_limit < 0:
        parser.error('--legacy-extension-limit must be nonnegative')
    try:
        prepare(args)
    except (ImportError, ValueError, RuntimeError, OSError, subprocess.CalledProcessError) as error:
        parser.exit(1, 'pyAmpliCol preparation failed: %s\n' % error)


if __name__ == '__main__':
    main()
