"""Test cli result module."""
from click.testing import CliRunner
from uwg.cli.simulate import simulate
from uwg.cli.validate import validate
from uwg.cli import viz

import os


def test_viz():
    runner = CliRunner()
    result = runner.invoke(viz)
    assert result.exit_code == 0
    assert result.output.startswith('vi')
    assert result.output.endswith('z!\n')


def test_model_validate():
    """Test uwg validation."""

    runner = CliRunner()
    input_uwg = './tests/json/uwg.json'
    result = runner.invoke(validate, ['model', input_uwg])
    assert result.exit_code == 0

    input_uwg = './tests/json/custom_uwg.json'
    result = runner.invoke(validate, ['model', input_uwg])
    assert result.exit_code == 0


def test_param_validate():
    """Test uwg validation."""

    runner = CliRunner()
    input_uwg = './tests/parameters/initialize_singapore.uwg'
    result = runner.invoke(validate, ['param', input_uwg])
    assert result.exit_code == 0

    input_uwg = './tests/parameters/initialize_singapore_error.uwg'
    result = runner.invoke(validate, ['param', input_uwg])
    assert result.exit_code == 1


def test_uwg_simulate():
    """Test simple uwg simulation."""

    runner = CliRunner()

    # test with basic command
    input_uwg = './tests/json/uwg.json'

    result = runner.invoke(simulate, ['model', input_uwg])
    assert result.exit_code == 0

    output_epw = './tests/epw/SGP_Singapore.486980_IWEC_UWG.epw'
    assert os.path.isfile(output_epw)
    os.remove(output_epw)

    # test with optional arguments
    new_epw_dir = './tests/epw_uwg'
    new_epw_name = 'custom_test_singapore.epw'
    epw_path = './tests/epw/SGP_Singapore.486980_IWEC.epw'
    result = runner.invoke(
        simulate, ['model', input_uwg, '--epw-path', epw_path, '--new-epw-dir',
                   new_epw_dir, '--new-epw-name', new_epw_name])
    assert result.exit_code == 0

    default_output_epw = './tests/epw/SGP_Singapore.486980_IWEC_UWG.epw'
    assert not os.path.isfile(default_output_epw)

    output_epw = os.path.join(new_epw_dir, new_epw_name)
    assert os.path.isfile(output_epw)
    os.remove(output_epw)


def test_custom_uwg_simulate():
    """Test custom uwg simulation."""

    runner = CliRunner()

    # test with basic command
    input_uwg = './tests/json/custom_uwg.json'
    epw_path = './tests/epw/SGP_Singapore.486980_IWEC.epw'
    result = runner.invoke(simulate, ['model', input_uwg, '--epw-path', epw_path])
    assert result.exit_code == 0

    output_epw = './tests/epw/SGP_Singapore.486980_IWEC_UWG.epw'
    assert os.path.isfile(output_epw)
    os.remove(output_epw)


def test_param_simulate():
    """Test simple uwg simulation from param file."""

    runner = CliRunner()

    # test with basic command
    input_uwg = './tests/parameters/initialize_singapore.uwg'
    epw_path = './tests/epw/SGP_Singapore.486980_IWEC.epw'

    result = runner.invoke(simulate, ['param', input_uwg, epw_path])
    assert result.exit_code == 0

    output_epw = './tests/epw/SGP_Singapore.486980_IWEC_UWG.epw'
    assert os.path.isfile(output_epw)
    os.remove(output_epw)
