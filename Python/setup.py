#!/usr/bin/env python3

import setuptools
from setuptools.extension import Extension
from Cython.Build import cythonize

clib = ("clib",
        {"sources": ["src/mini-gmp/mini-gmp.c"]
         })

ext_modules = [
    Extension("exact_tests",
              sources = ["src/exact_tests/exact.pyx",
                         "src/include/binom.cc",
                         "src/include/fisher.cc",
                         "src/include/plink2_base.cc",
                         "src/include/plink2_float.cc",
                         "src/include/plink2_hwe.cc",
                         "src/include/plink2_highprec.cc",
                         "src/include/plink2_ln.cc"],
              language = "c++",
              extra_compile_args = ["-std=c++11", "-Wno-unused-function", "-Wno-cpp", "-DNO_CPP11_TYPE_ENFORCEMENT"],
              extra_link_args = ["-std=c++11"]
              )
    ]

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="exact_tests",
    version="0.3.4",
    author="Christopher Chang",
    author_email="chrchang@alumni.caltech.edu",
    description="Accurate and efficient binomial, Hardy-Weinberg equilibrium, and Fisher's exact tests.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/chrchang/stats",
    project_urls={
        "Bug Tracker": "https://github.com/chrchang/stats/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU Lesser General Public License v3.0 only (LGPL-3.0-only)"
        "Operating System :: OS Independent",
    ],
    package_dir={"": "src"},
    packages=setuptools.find_namespace_packages(where="src"),
    python_requires=">=3.10",
    libraries=[clib],
    ext_modules=cythonize(ext_modules),
)
