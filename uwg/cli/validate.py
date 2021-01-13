
"""UWG JSON model validation commands."""

try:
    import click
except ImportError:
    raise ImportError(
        'click is not installed. Try `pip install . [cli]` command.'
    )

from uwg import UWG

import sys
import os
import json

try:
    import uwg_schema.model as schema_model
except ImportError:
    print('uwg_schema is not installed. Try `pip install . [cli]` command.')


@click.group(help='Commands for validating UWG JSON and .uwg files.')
def validate():
    pass


@validate.command('model')
@click.argument('model-json')
def validate_model(model_json):
    """Validate a UWG model JSON file against the UWG schema.

    \b
    Args:
        model_json: Full path to a UWG model JSON file. Note that this will
            not check if rural and new .epw file paths are valid since these
            can be overridden by CLI arguments.
    """
    try:
        assert os.path.isfile(
            model_json), 'No JSON file found at {}.'.format(model_json)

        # validate the Model JSON
        click.echo('Validating Model JSON ...')
        schema_model.UWG.parse_file(model_json)
        click.echo('Pydantic validation passed.')
        with open(model_json) as json_file:
            data = json.load(json_file)
        UWG.from_dict(data)
        click.echo('Python re-serialization passed.')
        click.echo('Congratulations! Your UWG model JSON is valid!')
    except Exception as e:
        print('Model validation failed.\n{}'.format(e))
        sys.exit(1)
    else:
        sys.exit(0)


@validate.command('param')
@click.argument('param-uwg')
def validate_param(param_uwg):
    """Validate a .uwg parameter file.

    \b
    Args:
        param_uwg: Full path to a .uwg parameter file.
    """
    try:
        # validate the Model JSON
        click.echo('Validating .uwg parameter file ...')
        UWG.from_param_file(param_uwg)
        click.echo('Congratulations! Your .uwg parameter file is valid!')
    except Exception as e:
        print('The .uwg parameter file validation failed.\n{}'.format(e))
        sys.exit(1)
    else:
        sys.exit(0)
