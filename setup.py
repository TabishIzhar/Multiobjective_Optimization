import pathlib
from setuptools import setup
from setuptools import find_packages

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
# README = (HERE / "README.md").read_text()

# This call to setup() does all the work
setup(
    name="Multiobjective optimization",
    version="0.1",
    description="",
    long_description= " ",
    long_description_content_type="text/markdown",
    url= " " ,
    author="Tabish Izhar",
    author_email="tizhar@iul.ac.in",
    license="MIT",
    classifiers=[
        "Programming Language :: Python :: 3.10.7",
    ],
    
    packages=find_packages(include=['MOSA', 'MOSA.*']),
    
    
    include_package_data=False,
    install_requires=["numpy", "pandas", "plotly"]
)