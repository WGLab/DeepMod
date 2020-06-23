#!/usr/bin/env python

import os,sys
import setuptools


with open("README.md", "r") as fh:
    long_description = fh.read()


setuptools.setup(
    name="NanoMod", # Replace with your own username
    version="0.1.2",
    author="Qian Liu",
    author_email="liuqianhn@gmail.com",
    description="A deep-learning tool to detect DNA modifications using Nanopore long-read sequencing data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/WGLab/DeepMod",
    #packages=setuptools.find_packages(),
    packages=['DeepMod_scripts'],
    package_dir={'DeepMod_scripts': 'bin/DeepMod_scripts'},
    package_data={'train_deepmod': ['train_deepmod/*']},
    classifiers=[
        "Programming Language :: Python",
        'Intended Audience :: Science/Research',
        "License :: OSI Approved :: GNU Library or Lesser General Public License (LGPL)",
        "Operating System :: OS Independent",
    ],
)
