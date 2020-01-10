#!/usr/bin/env python3

from distutils.core import setup

LONG_DESCRIPTION = \
'''
This script is designed to run a series of Unix bioinformatics tools in an
effort to assess the benefit of different alignment tools for metagenomic
binning, in particular assessment of contig depth of coverage.

'''


setup(
    name='MAGsim',
    version='0.1.0.0',
    author='Alistair Legione',
    author_email='legionea@unimelb.edu.au',
    packages=['MAGsim'],
    package_dir={'MAGsim': 'MAGsim'},
    entry_points={
        'console_scripts': ['MAGsim = MAGsim.pyContigDepthAnalysisPipeline:main']
    },
    url='https://github.com/alegione/MAGsim',
    license='MIT',
    description=('Simulation of metagenomic assemblies'),
    long_description=(LONG_DESCRIPTION),
    install_requires=["argparse", "os", "sys", "subprocess", "logging", "shutils",
    "pkg_resources", "pandas", "datetime", "Bio", "numpy"],
)
