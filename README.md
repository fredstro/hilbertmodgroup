# Hilbert Modular Groups

This repository contains a python package `hilbert_modgroup` that implements algorithms 
for Hilbert modular groups, in particular a reduction algorithm. The implementation is written in Python 
and is dependent on SageMath.

## Requirements
- SageMath v9.4 (https://www.sagemath.org/)
- Note: the main incompatibility with previous versions is that the 
  'Transformation' option to LLL was only introduced in v9.4. This means in particular that 
  the basic functionality likely work with 9.1-3 but the pre-optimization in finding a close cusp using LLL
  will fail. As a temporary workaround I have added a try/except to catch this.  

## Installation
With the executable `sage` in your path run 
   `make install`

## Usage
The package can be imported and used as any other package, e.g. 

```
sage: from hilbert_modgroup.all import HilbertModularGroup
sage: g = HilbertModularGroup(5) # Hilbert modular group for Q(sqrt(5))
```
For more examples see the embedded doctests (search for `EXAMPLES`) as well as
the `/examples` directory which contains Jupyter notebook with more extensive 
examples corresponding to the paper
"Reduction Algorithms for Hilbert Modular Groups" by F. Stromberg. 
(Reference to appear)

## Test
You can run all doctest by typing
`make test`

If you just want to test an individual module you can simply run e.g. 
`sage -t <file_to_test>`

## Docker
It is also possible to install the package and run the example notebooks
in a docker container with the following steps:
1. docker build -t hilbertmodgroup .
2. docker run -p 8888:8888 hilbertmodgroup
3. Open a browser at the indicated URL: http://127.0.0.1:8888/?token=<token>

