#!/usr/bin/env python3

from distutils.core import setup

LONG_DESCRIPTION = \
'''
This script is designed to run a series of Unix bioinformatics tools in an
effort to assess the benefit of different alignment tools for metagenomic
binning, in particular assessment of contig depth of coverage.

'''


setup(
    name='MAG-sim',
    version='0.1.0.0',
    author='Alistair Legione',
    author_email='legionea@unimelb.edu.au',
    packages=['MAG-sim'],
    package_dir={'MAG-sim': 'MAG-sim'},
    entry_points={
        'console_scripts': ['MAGsim = MAGsim.pyContigDepthAnalysisPipeline:main']
    },
    url='https://github.com/alegione/MAG-sim',
    license='MIT',
    description=('Simulation of metagenomic assemblies'),
    long_description=(LONG_DESCRIPTION),
    install_requires=["argparse", "os", "sys", "subprocess", "logging", "shutils",
    "pkg_resources", "pandas", "datetime", "Bio", "numpy"],
)
