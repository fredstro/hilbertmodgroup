# Hilbert Modular Groups

This repository contains a python package `hilbert_modgroup` that implements algorithms 
for Hilbert modular groups, in particular a reduction algorithm. The implementation is written in Python 
and is dependent on SageMath.

## Requirements
- SageMath v9.6+ (https://www.sagemath.org/)
  (Tested on v9.6, v9.7 and v10.0)

## Installation
### Using sage pip
This package needs to be installed in the virtual environment provided by SageMath, and it is therefore necessary 
to run the following command 
```console
$ sage -pip install --no-build-isolation hilbert-modular-group
```
**Note**: The `--no-build-isolation` is necessary as the compiler needs access 
to certain library files from the sage installation and SageMath itself is 
too large to be required as a build dependency. 
As an alternative to this flag you can also specify the environment variable 
SAGE_LIB explicitly.

### From git source
If the SageMath executable `sage` is in the current path you can install from source using the Makefile

```console
$ git clone https://github.com/fredstro/hilbertmodgroup.git
$ cd hilbertmodgrup
$ make install
```

### Docker
If you do not have SageMath installed, but you have docker you can use install this package
in a docker container built and executed using e.g. `make docker-sage` or `make docker-examples`


## Usage
The package can be imported and used as any other package. 
For example, to find the reduction of the point given by [1+i,1+i] in H^2 
with respect to the Hilbert modular group of Q joint by square-root of 5 write: 

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

## Examples

The directory `/examples` contains Jupyter notebooks with example code to illustrate the interface and functionality of this package. 
You can either open them manually from SageMath or run one of the following commands:
`make examples`
`make docker-examples`
which will start up a Jupyter notebook server from sagemath either locally or in a docker container. 

## Community Guidelines

### How to Contribute?
- Open an issue on GitHub and create a pull / merge request against the `develop` branch.
### How to report an issue or a problem? 
- First check if the issue is resolved in the `develop` branch. If not, open an issue on GitHub. 
### How to seek help and support?
- Contact the maintainer, Fredrik Stromberg, at: fredrik314@gmail.com (alternatively at fredrik.stromberg@nottingham.ac.uk)

## Development and testing

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
7. `tox` -- run `sage -tox` with all environments: `doctest`, `coverage`, `pycodestyle`, `relint`, `codespell`
   Note: If your local SageMath installation does not contain tox this will run `sage -pip install tox`.

The following commands are run in an isolated docker container 
and requires docker to be installed and running:
1. `docker` -- build a docker container with the tag `hilbertmodgroup-{GIT_BRANCH}`
2. `docker-rebuild` -- rebuild the docker container without cache
3. `docker-test` -- run SageMath's doctests in the docker container
4. `docker-examples` -- run a Jupyter notebook with the SageMath kernel initialised at the `/examples` directory 
  and exposing the notebook at http://127.0.0.1:8888. The port used can be modified by 
5. `docker-tox` -- run tox with all environments: `doctest`, `coverage`, `pycodestyle`, `relint`, `codespell`. 
6. `docker-shell` -- run a shell in a docker container
7. `docker-sage` -- run a sage interactive shell in a docker container

The following command-line parameters are available 
- `NBPORT` -- set the port of the notebook for `examples` and `docker-examples`  (default is 8888)
- `TOX_ARGS` -- can be used to select one or more of the tox environments (default is all)
- `REMOTE_SRC` -- set to 0 if you want to use the local source instead of pulling from gitHub (default 1)
- `GIT_BRANCH` -- the branch to pull from gitHub (used if REMOTE_SRC=1)

### Example usage
Run tox coverage on the branch `main` from gitHub:

`make docker-tox REMOTE_SRC=1 GIT_BRANCH=main TOX_ARGS=coverage`

Run doctests on the local source with local version of sage:

`make tox TOX_ARGS=doctest`

Run relint on the local source with docker version of sage:

`make docker-tox REMOTE_SRC=0 TOX_ARGS=relint`

## Development

### GitHub Workflow

- There are two long-lived branches `main` and `develop`.
- The `develop` branch is used for development and can contain new / experimental features.  
- Pull-requests should be based on `develop`.
- Releases should be based on `main`.
- The `main` branch should always be as stable and functional as possible. In particular, merges should always happen from `develop` into `main`. 
- Git-Flow is enabled (and encouraged) with feature branches based on `develop` and hotfixes based on `main`. 

### GitHub Actions

Each commit is tested and checked using gitHub actions with tox running:
- `doctest` -- run all doctests
- `coverage` -- ensure that all functions and classes are documented 
- `pycodestyle-minimal` -- ensure PEP8 style guide is followed (except we allow max line length 99)
- `relint` -- relint against some patterns taken from the SageMath source (config file .relint.yaml)
- `codespell` -- spellchecker

To make sure that your commit passes all tests you should `make tox` or `make docker-tox REMOTE_SRC=0` on the command line.

### Versions

Versioning of this project is managed by setuptools_scm.
To bump the version create a git tag `x.y.z` and the file 
 `src/hilbert_modgroup/version.py` will then be automatically updated to contain 
```
version = 'x.y.z.???'
version_tuple = (x, y, z, '???')
```
where ??? depends on the state of the current directory. 
If you are creating a new version to release the source directory should be clean.

### PyPi

To upload new versions to PyPi: 
1. `make sdist` -- creates a source distribution `dist/hilbert_modular_group-x.y.z`
2. `twine check dist/hilbert_modular_group-x.y.z`
3. `twine upload --repository pypi dist/hilbert_modular_group-z.y.z`

## References:

- [![DOI](https://joss.theoj.org/papers/10.21105/joss.03996/status.svg)](https://doi.org/10.21105/joss.03996)
