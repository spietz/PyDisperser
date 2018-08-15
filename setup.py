#!/usr/bin/env python
# distutils: language = c++
# distutils: sources = dispersion.cpp

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(
    ext_modules = [Extension(
        "mymod",   # name of extension
        ["mymod.pyx","dispersion.cpp"],  # filename of our Pyrex/Cython source
        language="c++",  # this causes Cython to create C++ source
        include_dirs=["m"],  # usual stuff
        libraries=["boost_python"],  # ditto
        extra_link_args=[],  # if needed
        extra_compile_args=["-std=c++0x"],
    )], cmdclass = {'build_ext': build_ext}
    
)
