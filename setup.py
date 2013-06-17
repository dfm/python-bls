#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
from numpy.distutils.core import setup, Extension


if sys.argv[-1] == "publish":
    os.system("python setup.py sdist upload")
    sys.exit()


# First, make sure that the f2py interfaces exist.
interface_exists = os.path.exists("bls/bls.pyf")
if "interface" in sys.argv or not interface_exists:
    # Generate the Fortran signature/interface.
    cmd = "cd bls;"
    cmd += "f2py eebls.f -m _bls -h bls.pyf"
    cmd += " --overwrite-signature"
    os.system(cmd)
    if "interface" in sys.argv:
        sys.exit(0)

# Define the Fortran extension.
bls = Extension("bls._bls", ["bls/bls.pyf", "bls/eebls.f"])

setup(
    name="bls",
    url="https://github.com/dfm/python-bls",
    version="0.0.1",
    author="Dan Foreman-Mackey",
    author_email="danfm@nyu.edu",
    description="",
    long_description="",
    packages=["bls", ],
    ext_modules=[bls, ],
    classifiers=[
        # "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
    ],
)
