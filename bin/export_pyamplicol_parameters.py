#!/usr/bin/env python
"""Read a MadNkLO model/card with its Python 2 interpreter for artifact preparation."""
import json
import os
import sys
import glob

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
from madgraph.iolibs import save_load_object
from models import import_ufo, model_reader, check_param_card


def export(directory):
    dump = save_load_object.load_from_file(os.path.join(directory, 'MadEvent7.db'))
    if dump['model_with_CMS']:
        raise ValueError('Complex-mass scheme is not supported by the provider adapter')
    return export_from_model(directory,dump['model_name'])


def export_from_model(directory, model_name):
    """Read a real exported model/card for bounded amplitude-only controls."""
    model_path = os.path.join(directory, 'Source', model_name)
    model = import_ufo.import_model(model_path, restrict=False, prefix=True, complex_mass_scheme=False)
    model = model_reader.ModelReader(model)
    card = check_param_card.ParamCard(os.path.join(directory,'Cards','param_card.dat'))
    defaults = None
    for path in glob.glob(os.path.join(directory,'*','BornProvider','*.json')):
        current = json.load(open(path)).get('restriction_defaults')
        if current is not None:
            if defaults is not None and defaults != current:
                raise ValueError('Conflicting model restrictions in provider manifests')
            defaults = current
    if defaults is None:
        raise ValueError('A fresh export is required to preserve the original model restriction')
    for block,code,value in defaults:
        if not card.has_param(block,code):
            card.add_param(block,code,value)
    model.set_parameters_and_couplings(
        param_card=card,
        complex_mass_scheme=False)
    parameters = {}
    for parameter in model.get('parameters')[('external',)]:
        value = complex(model.get('parameter_dict')[parameter.name])
        name = parameter.name[4:] if parameter.name.startswith('mdl_') else parameter.name
        parameters[name] = [value.real, value.imag]
    return {'model': model_path, 'parameters': parameters,'restriction_card':card.write()}


if __name__ == '__main__':
    # Imported libraries can log to stdout; the marker delimits the JSON payload.
    print('PYAMP_METADATA:' + json.dumps(export(os.path.abspath(sys.argv[1]))))
