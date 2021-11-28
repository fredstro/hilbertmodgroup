
import os
import setuptools

from setuptools.extension import Extension
from Cython.Build import cythonize

from sage.env import sage_include_directories,SAGE_LIB,SAGE_LOCAL,NTL_INCDIR,NTL_LIBDIR
SAGE_INC=SAGE_LOCAL + "/include"

include_path=[SAGE_INC,SAGE_LIB]

# Extension modules using Cython
ext_modules=[
    Extension(
        'hilbert_modgroup.hilbert_modular_group_element',
        sources=[os.path.join('src/hilbert_modgroup/hilbert_modular_group_element.pyx')]
    ),
    Extension(
        'hilbert_modgroup.upper_half_plane',
        sources=[os.path.join('src/hilbert_modgroup/upper_half_plane.pyx')],
        extra_compile_args=['-Wno-unused-function']
    ),
    Extension(
        'hilbert_modgroup.pullback_cython',
        sources=[os.path.join('src/hilbert_modgroup/pullback_cython.pyx')],
        include_dirs=[os.path.join(SAGE_LIB,'cypari2')],
        language='c++',
        extra_compile_args=['-std=c++11'])
    ,
]

setuptools.setup(
    name='hilbert_modular_group',
    version='1.0',
    license='GPL v3+',
    author='Fredrik Stromberg',
    author_email='fredrik314@gmail.com',
    packages=['hilbert_modgroup'],
    package_dir={'hilbert_modgroup':'src/hilbert_modgroup'},
    zip_safe=False,
    ext_modules=cythonize(
        ext_modules,
        include_path=['src'],
        compiler_directives={
            'embedsignature': True,
            'language_level': '3',
        },
    ),
)
