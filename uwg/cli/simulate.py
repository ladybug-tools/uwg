"""uwg simulation running commands."""

try:
    import click
except ImportError:
    raise ImportError(
        'click is not installed. Try `pip install . [cli]` command.'
    )

from uwg import UWG

import sys
import json


@click.group(help='Commands for simulating UWG models.')
def simulate():
    pass


@simulate.command('model')
@click.argument('model-json', type=click.Path(
    exists=True, file_okay=True, dir_okay=False, resolve_path=True))
@click.argument('epw-path', type=click.Path(
    exists=True, file_okay=True, dir_okay=False, resolve_path=True))
@click.option('--new-epw-dir', help='Optional argument for the destination directory '
              'into which the morphed .epw file is written. The argument passed here '
              'will overwrite the new_epw_dir specified in the UWG JSON model file.',
              default=None, show_default=True)
@click.option('--new-epw-name', help='Optional argument for The destination file name '
              'of the morphed .epw file. The argument passed here will overwrite the '
              'new_epw_name specified in the UWG JSON model file.',
              default=None, show_default=True)
def simulate_json_model(model_json, epw_path, new_epw_dir, new_epw_name):
    """Simulate a UWG model from a JSON model file.

    \b
    Args:
        model_json: Full path to a JSON model file.
        epw_path: Full path of the rural .epw file that will be morphed.
    """
    try:
        with open(model_json) as json_file:
            data = json.load(json_file)

        uwg_model = UWG.from_dict(
            data, epw_path=epw_path, new_epw_dir=new_epw_dir, new_epw_name=new_epw_name)
        uwg_model.generate()
        uwg_model.simulate()
        uwg_model.write_epw()
    except Exception as e:
        print('UWG model simulation failed.\n{}'.format(e))
        sys.exit(1)
    else:
        sys.exit(0)


@simulate.command('param')
@click.argument('param-uwg', type=click.Path(
    exists=True, file_okay=True, dir_okay=False, resolve_path=True))
@click.argument('epw-path',  type=click.Path(
    exists=True, file_okay=True, dir_okay=False, resolve_path=True))
@click.option('--new-epw-dir', help='Optional argument for the destination directory '
              'into which the morphed .epw file is written.',
              default=None, show_default=True)
@click.option('--new-epw-name', help='Optional argument for The destination file name '
              'of the morphed .epw file.',
              default=None, show_default=True)
def simulate_uwg_param_model(param_uwg, epw_path, new_epw_dir, new_epw_name):
    """Simulate a UWG model from a .uwg parameter file.

    \b
    Args:
        param_uwg: Full path to a .uwg param file.
        epw_path: Full path of the rural .epw file that will be morphed.
    """
    try:
        uwg_model = UWG.from_param_file(
            param_uwg, epw_path, new_epw_dir, new_epw_name)
        uwg_model.generate()
        uwg_model.simulate()
        uwg_model.write_epw()
    except Exception as e:
        print('UWG model simulation failed.\n{}'.format(e))
        sys.exit(1)
    else:
        sys.exit(0)
