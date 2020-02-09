# -*- coding: utf-8 -*-

# Learn more: https://github.com/kennethreitz/setup.py

from setuptools import setup, find_packages


with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='sc_topics',
    version='0.1.0',
    description='Interface to topic modeling for single cell RNA-seq',
    long_description=readme,
    author='Kirk Gosik',
    author_email='kgosik@broadinstitute.org',
    url='https://github.com/broadinstitute/sc_topics',
    license=license,
    packages=find_packages(exclude=('tests', 'docs'))
)