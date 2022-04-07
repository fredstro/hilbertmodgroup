import os
import setuptools
from setuptools.extension import Extension
from Cython.Build import cythonize
# Check if we are currently in a SageMath environment.
SAGE_LOCAL = os.getenv('SAGE_LOCAL')
if not SAGE_LOCAL:
    raise ValueError("This package can only be installed inside SageMath (http://www.sagemath.org)")
# Find correct value for SAGE_LIB which is needed to compile the Cython extensions.
SAGE_LIB = os.getenv('SAGE_LIB')
if not SAGE_LIB:
    try:
        from sage.env import SAGE_LIB
    except ModuleNotFoundError:
        raise ModuleNotFoundError("To install this package you need to either specify the "
                                  "environment variable 'SAGE_LIB' or call pip with "
                                  "'--no-build-isolation'")
if not os.path.isdir(SAGE_LIB):
    raise ValueError(f"The library path {SAGE_LIB} is not a directory.")

# Extension modules using Cython
extra_compile_args = ['-Wno-unused-function',
                      '-Wno-implicit-function-declaration',
                      '-Wno-unused-variable',
                      '-Wno-deprecated-declarations',
                      '-Wno-deprecated-register']
ext_modules = [
    Extension(
        'hilbert_modgroup.hilbert_modular_group_element',
        sources=[os.path.join('src/hilbert_modgroup/hilbert_modular_group_element.pyx')],
        extra_compile_args=extra_compile_args
    ),
    Extension(
        'hilbert_modgroup.upper_half_plane',
        sources=[os.path.join('src/hilbert_modgroup/upper_half_plane.pyx')],
        extra_compile_args=extra_compile_args
    ),
    Extension(
        'hilbert_modgroup.pullback_cython',
        sources=[os.path.join('src/hilbert_modgroup/pullback_cython.pyx')],
        language='c++',
        extra_compile_args=extra_compile_args+['-std=c++11'])
]

setuptools.setup(
    ext_modules=cythonize(
        ext_modules,
        include_path=['src',SAGE_LIB],
        compiler_directives={
            'embedsignature': True,
            'language_level': '3',
        },
    ),
)
