import setuptools
import os


setuptools.setup(
    name='betazoid',
    version='0.1',    
    description='N/A',
    author='Philippa Richter',
    author_email='prichter@berkeley.edu',
    packages=['src', 'src.files'],
    entry_points={'console_scripts': ['recruit_reads=src.cli:recruit_reads', 'align_reads=src.cli:recruit_reads']})

# TODO: What exactly is an entry point?
# https://python-packaging.readthedocs.io/en/latest/command-line-scripts.html
#  