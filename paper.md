---
title: 'hilbertmodgroup: Reduction algorithms and framework for Hilbert Modular Groups' 

tags:
  - Python
  - mathematics
  - SageMath
authors:
  - name: Fredrik Strömberg 
    orcid:0000-0002-6960-1002
    affiliation: "University of Nottingham"

date: 29 November 2021
bibliography: paper.bib

---

# Summary

This package implements basic classes and a new reduction algorithm for Hilbert modular groups. 
The main improvement over previous algorithms is that
this implementation works in theory for all Hilbert modular groups and in practice for a much 
wider range of examples. A more in-depth discussion of the theoretical background and details about the implementation can be found in [@reductionalghilbert].


# A Brief Mathematical Background 

One of the most important groups in Number Theory is the modular group, $\Gamma=\mathrm{PSL}_2(\mathbb{Z})$, 
consisting of fractional-linear transformations, $z\mapsto (az+b)/(cz+d)$ on the complex upper
half-plane, $\mathbb{H}$, given by 2-by-2 matrices of determinant 1 and integer entries. 

A reduction algorithm theory for the modular group is an algorithm that, for a given $z \in \mathbb{H}$, 
finds an element, $A=\begin{smallmatrix} a & b \\ c & d\end{smallmatrix} \in \Gamma$ such that $Az=(az+b)/(cz+d)$ 
belongs to a specific set, a so-called "fundamental domain". This type of algorithm was first 
found by Gauss in connection with binary quadratic forms.

A natural generalisation of the modular group over $\mathbb{Z}$, 
is given by the family of Hilbert modular groups, $\Gamma_K=\mathrm{PSL}_2(\mathcal{O}_K)$, 
where $K$ is a totally real number field of degree $n$ and $\mathcal{O}_K$ is its ring of integers. 
This group gives rise to an action on $n$ copies of the complex upper half-plane 
$$\mathbb{H}_K=\mathbb{H} \times \cdots \times \mathbb{H}.$$ 
A reduction algorithm for a Hilbert modular group $\Gamma_K$ 
works in the same way as before. Given $z \in \mathbb{H}_K$ the algorithm finds an element 
$A \in \Gamma_K$ such that $Az$ belongs to a certain fundamental domain. 
The additional complexity it in this case, when $K$ is not equal to $\QQ$,
has both a theoretical and a numerical part. stems both from 

# Statement of need

There has been several previous attempts at giving a reduction algorithm for Hilbert modular groups 
but they have all been limited in at least one of two ways with 
 the number field either being restricted to degree 2, or the class number to be $1$, or a combination of both. 
See for example the algorithms by Boyer and Streng [@MR3376741], and Quinn and Verjovsky [@MR4091535].


Having access to the algorithm in this package, which is valid for any totally real number field, 
opens up for several new research directions and generalisations of previous research. 
Some of the direct applications lie in the field of explicit formulas for 
Hilbert modular groups and computational aspects of non-holomorphic Hilbert modular forms. 

# Implementation

The package `hilbertmodgrup` is mainly written in Python with some parts in Cython [@behnel2011cython]. 
It is intended to run as a package inside SageMath [@sage] as it 
makes heavy use SageMath's implementation of 
number fields, which is in turn is in many cases using the backend from PARI/gp [@PARI2].

# Documentation and Testing
All functions are documented using docstrings with integrated doctest 
following the guide for SageMath development. 
In addition, the directory `/examples` contains a selection of 
Jupyter notebooks illustrating the use of the package with 
a selection of fundamental examples, corresponding to examples listed in 
[@reductionalghilbert].


# Acknowledgements
This work was partially supported by the Engineering and Physical Sciences Research Council [EP/N007360/1] and [EP/V026321/1].

# References