# Hilbert Modular Groups

This repository contains a python package `hilbert_modgroup` that implements algorithms 
for Hilbert modular groups, in particular a reduction algorithm. The implementation is written in Python 
and is dependent on SageMath.

## Requirements
- SageMath v9.4+ (https://www.sagemath.org/)

## Installation
### Using sage pip
This package needs to be installed in the virtual environment provided by SageMath and it is therefore necessary 
to run the install command 
```console
$ sage -pip install hilbert-modular-group
```

### From git source
If the SageMath executable `sage` is in the current path you can install from source using the Makefile

```console
$ git clone https://github.com/fredstro/hilbertmodgroup.git
$ cd hilbertmodgrup
$ make install
```

## Usage
The package can be imported and used as any other package. 
For example, to find the reduction of the point given by [1+i,1+i] in H^2 
with respect to the Hilbert modular group of Q adjoint square-root of 5 write: 

```
sage: from hilbert_modgroup.all import *
sage: H1=HilbertModularGroup(5)
sage: P1=HilbertPullback(H1)
sage: z = UpperHalfPlaneProductElement([1+I,1+I])
sage: P1.reduce(z)
[1.00000000000000*I, 1.00000000000000*I]
sage: z = UpperHalfPlaneProductElement([0.25+I/2,1+I])
sage: P1.reduce(z) # abs tol 1e-10
[0.694427190999916 + 0.611145618000168*I, -0.309016994374947 + 1.30901699437495*I]
sage: P1.reduce(z, return_map=True)[1]
[-1/2*a + 1/2  1/2*a + 1/2]
[-1/2*a + 1/2            0]

```
For more examples see the embedded doctests (search for `EXAMPLES`) as well as
the `/examples` directory which contains Jupyter notebook with more extensive 
examples corresponding to the paper
"Reduction Algorithms for Hilbert Modular Groups" by F. Stromberg. (Reference to appear)

## Additional Commands
The make file `Makefile` contains a number of useful commands that you can run using 
```console
$ make <command>
```
The following commands are run in your local SagMath environment:
1. `build` -- builds the package in place (sometimes useful for development).
2. `sdist` -- create a source distribution in /sdist (can be installed using `sage -pip install sdist/<dist name>`)
3. `install` -- build and install the package in the currently active sage environment
4. `clean` -- remove all build and temporary files
5. `test` -- run sage's doctests (same as `sage -t src/*`)
6. `examples` -- run a Jupyter notebook with the SageMath kernel initialised at the `/examples` directory.
7. `tox` -- run tox with all environments: `doctest`, `coverage`, `pycodestyle-minimal`, `relint`, `codespell`

The following commands are run in an isolated docker container 
and requires docker to be installed and running:
1. `docker` -- build a docker container with the tag `hilbertmodgroup`
2. `docker-rebuild` -- rebuild the docker container without cache
3. `docker-test` -- run SageMath's doctests in the docker container
4. `docker-examples` -- run a Jupyter notebook with the SageMath kernel initialised at the `/examples` directory 
  and exposing the notebook at http://127.0.0.1:8888
5. `docker-tox` -- run tox with all environments: `doctest`, `coverage`, `pycodestyle-minimal`, `relint`, `codespell`
6. `docker-shell` -- run a shell in a docker container
7. `docker-sage` -- run a sage interactive shell in a docker container

## Development

Each commit is tested and checked using github actions with tox running:
- `doctest` -- run all doctests
- `coverage` -- ensure that all functions and classes are documented 
- `pycodestyle-minimal` -- ensure PEP8 style guide is followed (except we allow max line length 99)
- `relint` -- relint against some of the patterns taken from the SageMath source (config file .relint.yaml)
- `codespell` -- spellchecker

To make sure that your commit passes all tests you can run `make tox` or `make docker-tox` on the command line. 