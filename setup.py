#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 Malcolm Ramsay <malramsay64@gmail.com>
#
# Distributed under terms of the MIT license.

# Always prefer setuptools over distutils
from setuptools import find_packages, setup

setup(
    name="dynamics_analysis",
    version="0.1.0",
    description="Analysis of dynamics",
    url="https://github.com/malramsay64/Dynamics",
    license="LICENSE",
    author="Malcolm Ramsay",
    author_email="malramsay64@gmail.com",  # Optional
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
    ],
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    python_requires=">=3.6,<4",
    install_requires=[],
    extras_require={"dev": [], "test": [],},
)
