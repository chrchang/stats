#!/usr/bin/env python3

import setuptools
from setuptools.extension import Extension
from Cython.Build import cythonize

# Note that this package makes heavy use of high-precision arithmetic functions
# which get a huge speedup from Fused Multiply-Add (FMA) instructions,
# available in x86-64-v3 but not base x86-64.  As a consequence, once PEP 817
# is finalized (and possibly even before that), it is critical to take
# advantage of it.  If other build backends add support for PEP 817 before
# setuptools does, we must switch.

ext_modules = [
    Extension("exact_tests",
              sources = ["src/exact_tests/exact.pyx",
                         "src/include/binom.cc",
                         "src/include/binom_detail.cc",
                         "src/include/fisher.cc",
                         "src/include/hypergeom.cc",
                         "src/include/hypergeom_detail.cc",
                         "src/include/nchypergeom_fisher.cc",
                         "src/include/plink2_base.cc",
                         "src/include/plink2_float.cc",
                         "src/include/plink2_hwe.cc",
                         "src/include/plink2_highprec.cc",
                         "src/include/special_func.cc"],
              language = "c++",
              extra_compile_args = ["-std=c++14", "-Wno-unused-function", "-Wno-cpp", "-DNO_CPP11_TYPE_ENFORCEMENT", "-ffp-contract=off"],
              extra_link_args = ["-std=c++14"]
              )
    ]

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="exact_tests",
    version="0.7.1",
    author="Christopher Chang",
    author_email="chrchang@alumni.caltech.edu",
    description="Accurate and efficient binomial, Hardy-Weinberg equilibrium, and Fisher's exact tests, along with associated distributions.",
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
    ext_modules=cythonize(ext_modules),
)
