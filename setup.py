#!/bin/env python3

from setuptools import setup, find_packages
setup(
    name="HatsPil",
    version="0.3rc0",
    packages=find_packages(),
    install_requires=[
        "cutadapt",
        "formatizer",
        "rpy2",
        "pandas",
        "tables",
        "numpy",
        "PyVCF"
    ],
    extras_require={
        "MongoDB": ["pymongo"]
    },
    package_data={'': ["data.hdf"]},
    include_package_data=True,
    entry_points={
        "console_scripts": [
            "hatspil = hatspil.hatspil:main",
        ],
    },
    author="Edoardo Morandi",
    author_email="morandidodo@gmail.com",
    description="HTS pipeline for different purposes",
    license="GPL3"
)
