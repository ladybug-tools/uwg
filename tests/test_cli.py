"""Test cli result module."""
from click.testing import CliRunner
from uwg.cli.simulate import simulate
from uwg.cli.validate import validate
from uwg.cli import viz

import os
import sys


def test_viz():
    runner = CliRunner()
    result = runner.invoke(viz)
    assert result.exit_code == 0, result.output
    assert result.output.startswith('vi'), result.output
    assert result.output.endswith('z!\n'), result.output


def test_model_validate():
    """Test uwg validation."""

    if (sys.version_info >= (3, 7)):
        runner = CliRunner()
        input_uwg = './tests/json/uwg.json'
        result = runner.invoke(validate, ['model', input_uwg])
        assert result.exit_code == 0, result.output

        input_uwg = './tests/json/custom_uwg.json'
        result = runner.invoke(validate, ['model', input_uwg])
        assert result.exit_code == 0, result.output


def test_param_validate():
    """Test uwg validation."""

    runner = CliRunner()
    input_uwg = './tests/parameters/initialize_singapore.uwg'
    result = runner.invoke(validate, ['param', input_uwg])
    print(result.output)
    assert result.exit_code == 0, result.output

    input_uwg = './tests/parameters/initialize_singapore_error.uwg'
    result = runner.invoke(validate, ['param', input_uwg])
    assert result.exit_code == 1, result.output


def test_uwg_simulate():
    """Test simple uwg simulation."""

    runner = CliRunner()

    # test with basic command
    input_uwg = './tests/json/uwg.json'
    epw_path = './tests/epw/SGP_Singapore.486980_IWEC.epw'

    result = runner.invoke(simulate, ['model', input_uwg, epw_path])
    assert result.exit_code == 0, result.output

    output_epw = './tests/epw/SGP_Singapore.486980_IWEC_UWG.epw'
    assert os.path.isfile(output_epw)
    os.remove(output_epw)

    # test with optional arguments
    new_epw_dir = './tests/epw_uwg'
    new_epw_name = 'custom_test_singapore.epw'
    result = runner.invoke(
        simulate, ['model', input_uwg, epw_path, '--new-epw-dir',
                   new_epw_dir, '--new-epw-name', new_epw_name])
    assert result.exit_code == 0, result.output

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
    result = runner.invoke(
        simulate, ['model', input_uwg, epw_path])
    assert result.exit_code == 0, result.output

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
    assert result.exit_code == 0, result.output

    output_epw = './tests/epw/SGP_Singapore.486980_IWEC_UWG.epw'
    assert os.path.isfile(output_epw)
    os.remove(output_epw)
