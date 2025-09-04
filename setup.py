import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

with open('cli-requirements.txt') as f:
    cli_requirements = f.read().splitlines()

setuptools.setup(
    name="uwg",
    use_scm_version=True,
    setup_requires=['setuptools_scm'],
    author="Ladybug Tools",
    author_email="info@ladybug.tools",
    description="Python application for modeling the urban heat island effect.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ladybug-tools/uwg",
    packages=setuptools.find_packages(exclude=["tests*", "resources*"]),
    include_package_data=True,
    install_requires=[],
    extras_require={
        'cli': cli_requirements
    },
    entry_points={
        "console_scripts": ["uwg = uwg.cli:main"]
    },
    classifiers=[
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: Implementation :: CPython",
        "Operating System :: OS Independent"
    ],
    license="MIT"
)
