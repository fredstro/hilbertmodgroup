# Hilbert Modular Forms

This repository contains a python package `hilbert_modgroup` that implements algorithms 
for Hilbert modular forms. The implementation is written in Python 
and is dependent on SageMath.

## Requirements
- SageMath (https://www.sagemath.org/)
  Tested on v 9.4 but should work on most 9.+ versions.

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