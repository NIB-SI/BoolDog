#!/usr/bin/env python

from setuptools import setup

setup(
    name='BoolDoG',
    version="0.0.1",
    packages=['booldog'],
    python_requires=">=3.6",
    install_requires=[
        'pyboolnet @ git+https://github.com/hklarner/pyboolNet@3.0.1',
        'xmltodict>=0.12',
        'matplotlib>=3.4',
        'numpy>=1.20',
        'scipy>=1.6',
        'igraph>=0.9.11',
        'biosimulators-utils[logging,sbml]',
        'xmltodict>=0.13.0'
    ],
    extras_require={
        "SBML": ["python-libsbml>=5.19"],
        "networks": ["networkx", "igraph"],
        "graphviz": ["pygraphviz"],
        "biomodels": ["biomodels-restful-api-client>=0.1.3"]
    },
    author='Carissa Bleker',
    author_email='carissa.bleker@nib.si',
    url='https://github.com/NIB-SI/BoolDoG',
    description=
    'A Python package for analyses of Boolean and semi-quantitative Boolean networks.',
    license='GPL3',
)
