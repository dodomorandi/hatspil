#!/bin/env python3

from setuptools import find_packages, setup

setup(
    name="HatsPil",
    version="0.6",
    packages=find_packages(),
    python_requires=">=3.6",
    install_requires=[
        "cutadapt",
        "rpy2",
        "pandas",
        "tables",
        "numpy",
        "PyVCF",
        "plotly",
        "colorlover",
        "spectra",
    ],
    extras_require={"MongoDB": ["pymongo"]},
    package_data={"": ["data.hdf"]},
    include_package_data=True,
    entry_points={"console_scripts": ["hatspil = hatspil.hatspil:main"]},
    author="Edoardo Morandi",
    author_email="morandidodo@gmail.com",
    description="HTS pipeline for different purposes",
    license="MIT",
)
