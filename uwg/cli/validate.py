
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
import logging
import json

_logger = logging.getLogger(__name__)

try:
    import uwg_schema.model as schema_model
except ImportError:
    _logger.exception(
        'uwg_schema is not installed. Try `pip install . [cli]` command.'
    )


@click.group(help='Commands for validating UWG JSON files.')
def validate():
    pass


@validate.command('model')
@click.argument('model-json')
def validate_model(model_json):
    """Validate a UWG model JSON file against the UWG schema.
    \n
    Args:
        model_json: Full path to a UWG model JSON file.
    """
    try:
        assert os.path.isfile(model_json), 'No JSON file found at {}.'.format(model_json)

        # validate the Model JSON
        click.echo('Validating Model JSON ...')
        schema_model.UWG.parse_file(model_json)
        click.echo('Pydantic validation passed.')
        with open(model_json) as json_file:
            data = json.load(json_file)
        # overwrite epw_path so that it's not checked for validity.
        # since users may choose to override the path with cli
        data['epw_path'] = '.'
        UWG.from_dict(data)
        click.echo('Python re-serialization passed.')
        click.echo('Congratulations! Your UWG model JSON is valid!')
    except Exception as e:
        _logger.exception('Model validation failed.\n{}'.format(e))
        sys.exit(1)
    else:
        sys.exit(0)
