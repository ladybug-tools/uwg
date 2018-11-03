import re
import setuptools
import sys

with open("README.md", "r") as fh:
    long_description = fh.read()


with open("README.md", "r") as fh:
    long_description = fh.read()

with open('uwg/__init__.py', 'r') as fd:
    version = re.search(
        r'^__version__\s*=\s*[\'"]([^\'"]*)[\'"]',
        fd.read(),
        re.MULTILINE
    ).group(1)

try:
    from semantic_release import setup_hook
    setup_hook(sys.argv)
except ImportError:
    pass

setuptools.setup(
    name="uwg",
    version=version,
    author="ladybug-tools",
    description="The Urban Weather Generator engine for Urban Heat Island modelling",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ladybug-tools/uwg",
    packages=['uwg'],
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3.6",
        "Operating System :: OS Independent",
    ],
)
