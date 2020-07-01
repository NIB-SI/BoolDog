#!/usr/bin/env python

from setuptools import setup
import os


setup(
    name='squad-reboot',
    version="0.0.1",
    packages=['squad_reboot'],
    install_requires=['PyBoolNet', 'xmltodict', 'matplotlib', 'numpy', 'scipy', 'python-igraph'],

    author='Carissa Bleker',
    author_email='carissarobyn.bleker@nib.si',
    url='https://github.com/NIB-SI/squad-reboot', 
    description='Python implemetation of SQUAD algorithm', 
    license='MIT'
)