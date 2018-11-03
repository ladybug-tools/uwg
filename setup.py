import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="uwg",
    version="5.0.0",
    author="ladybug-tools",
    description="The Urban Weather Generator engine for Urban Heat Island modelling",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ladybug-tools/uwg",
    packages=setuptools.find_packages(),
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3.6",
        "Operating System :: OS Independent",
    ],
)
