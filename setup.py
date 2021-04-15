#!/usr/bin/env python

from setuptools import setup
import os


setup(
    name='squad-reboot',
    version="0.0.1",
    packages=['squad_reboot'],
    install_requires=[
    	'PyBoolNet @ git+https://github.com/hklarner/PyBoolNet@2.31.0',
    	'xmltodict>=0.12',
    	'matplotlib>=3.4',
    	'numpy>=1.20',
    	'scipy>=1.6',
    	'python-igraph>=0.9.1'
    	],
    author='Carissa Bleker',
    author_email='carissarobyn.bleker@nib.si',
    url='https://github.com/NIB-SI/squad-reboot',
    description='Python implemetation of SQUAD algorithm',
    license='MIT'
)