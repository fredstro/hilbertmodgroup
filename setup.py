import os
import sys
import setuptools
from setuptools.extension import Extension
from Cython.Build import cythonize

# Find correct value for SAGE_LIB without importing sage (to allow for build isolation)
SAGE_LOCAL = os.getenv('SAGE_LOCAL')
if not SAGE_LOCAL:
    raise ValueError("This package can only be installed inside SageMath (http://www.sagemath.org)")
# Need to find the correct site-packages directory
PYTHON_VERSION = f"{sys.version_info.major}.{sys.version_info.minor}"
SAGE_LIB = f"{SAGE_LOCAL}/lib/python{PYTHON_VERSION}/site-packages"
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
        include_path=[SAGE_LIB,'src'],
        compiler_directives={
            'embedsignature': True,
            'language_level': '3',
        },
    ),
)
