#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

with open("README.rst") as readme_file:
    readme = readme_file.read()

with open("HISTORY.rst") as history_file:
    history = history_file.read()

requirements = [
    "numpy",
    "scipy",
    "matplotlib",
    "ipython",
    "tqdm",
    "svgwrite",
    "runman",
    "ipython",
]

setup_requirements = ["pytest-runner"]

test_requirements = [
    "pytest",
    "pytest-cov",
    "pytest-raises",
    "codecov",
    "flake8",
    "black",
]

dev_requirements = [
    "black",
    "bumpversion>=0.5.3",
    "tox>=3.5.2",
    "coverage>=5.0a4",
    "Sphinx>=2.0.0b1",
    "sphinx_rtd_theme>=0.4.0",
    "twine>=1.13.0",
    "pytest>=4.3.0",
    "pytest-cov==2.6.1",
    "pytest-raises>=0.10",
    "pytest-runner>=4.4",
]

extra_requirements = {
    "test": test_requirements,
    "setup": setup_requirements,
    "dev": dev_requirements,
    "all": [*requirements, *test_requirements, *setup_requirements, *dev_requirements],
}

setup(
    author="Dave Williams",
    author_email="cdave@uw.edu",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Developers",
        "License :: Free for non-commercial use",
        "Natural Language :: English",
        "Programming Language :: Python :: 3.7",
    ],
    description="  ",
    install_requires=requirements,
    license="Allen Institute Software License",
    long_description=readme + "\n\n" + history,
    include_package_data=True,
    keywords="flins",
    name="flins",
    packages=find_packages(exclude=["tests", "*.tests", "*.tests.*"]),
    setup_requires=setup_requirements,
    test_suite="tests",
    tests_require=test_requirements,
    extras_require=extra_requirements,
    url="https://github.com/AllenCellModeling/flins",
    version="0.1.0",
    zip_safe=False,
)
