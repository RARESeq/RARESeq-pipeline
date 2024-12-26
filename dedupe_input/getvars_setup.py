# Setup file for getvars. This is a heavy function that is cythonized (cython) for
# performance. It takes in getvars.pyx and generates getvars.c and getvars.so. It 
# can be run using:
# python getvars_setup.py build_ext --inplace
from distutils.core import setup
from Cython.Build import cythonize

setup(
    ext_modules = cythonize('getvars.pyx')
)