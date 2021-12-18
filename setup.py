#!/usr/bin/env python

import os,sys
import setuptools


with open("README.md", "r") as fh:
    long_description = fh.read()


setuptools.setup(
    name="DeepMod", # Replace with your own username
    version="0.1.3",
    author="Qian Liu",
    author_email="",
    description="A deep-learning tool to detect DNA modifications using Nanopore long-read sequencing data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/WGLab/DeepMod",
    #packages=setuptools.find_packages(),
    packages=['DeepMod_scripts'],
    package_dir={'DeepMod_scripts': 'bin/DeepMod_scripts'},
    scripts=['bin/DeepMod.py', 'DeepMod_tools/cal_EcoliDetPerf.py', 'DeepMod_tools/generate_motif_pos.py', 'DeepMod_tools/hm_cluster_predict.py', 'DeepMod_tools/sum_chr_mod.py'],
    package_data={'train_deepmod': ['train_deepmod/*/*']},
    #data_files=[('train_deepmod', ['train_deepmod/*'])],
    classifiers=[
        "Programming Language :: Python",
        'Intended Audience :: Science/Research',
        "License :: OSI Approved :: GNU Library or Lesser General Public License (LGPL)",
        "Operating System :: OS Independent",
    ],
)
