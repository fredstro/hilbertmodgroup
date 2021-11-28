---
title: 'hilbertmodgroup: Reduction algorithms and framework for Hilbert Modular Groups'
tags:
  - Python
  - SageMath
  - mathematics
  - number theory
authors:
  - name: Fredrik Strömberg 
    orcid: 0000-0002-6960-1002
    affiliation: 1
affiliations:
  - name: "University of Nottingham, School of Mathematics"
    index: 1
date: 28 November 2021
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

A *reduction* algorithm for the modular group is an algorithm that, for a given $z \in \mathbb{H}$, 
finds an element, $A=\left(\begin{smallmatrix} a & b \\ c & d\end{smallmatrix}\right) \in \Gamma$ such that $Az=(az+b)/(cz+d)$ 
belongs to a specific set, a so-called "fundamental domain". 
This type of algorithm was first introduced in the context of binary quadratic forms in the 
18th century by Lagrange, Gauss and others, with the main contribution 
published by Gauss in his famous *Disquisitiones Arithmeticae*.

A natural generalisation of the modular group over $\mathbb{Z}$
is given by the family of Hilbert modular groups, $\Gamma_K=\mathrm{PSL}_2(\mathcal{O}_K)$, 
where $K$ is a totally real number field of degree $n$ and $\mathcal{O}_K$ is its ring of integers. 
This group gives rise to an action on $n$ copies of the complex upper half-plane 
$$\mathbb{H}_K=\mathbb{H} \times \cdots \times \mathbb{H}.$$ 
A reduction algorithm for a Hilbert modular group $\Gamma_K$ 
should work in the same way as before. Given $z \in \mathbb{H}_K$ the algorithm finds an element 
$A \in \Gamma_K$ such that $Az$ belongs to a certain fundamental domain. 
The additional complexity in this case, when $K$ is not equal to $\mathbb{Q}$,
has both a theoretical and a practical part. The main theoretical problem arises 
when the number field $K$ has class number greater than 1, in which case the corresponding
Fundamental domain will have more than one point at "infinity". 
From a practical standpoint the main problem appears when the degree and discriminant 
of the number field increases, making it necessary to, for instance, locate 
integral points in higher-dimensional polytopes. 

# Statement of need

There has been several previous attempts at giving a reduction algorithm for Hilbert modular groups 
but they have all been limited in at least one of two ways with 
 the number field either being restricted to degree 2, or the class number to be $1$, or a combination of both. 
See for example the algorithms by Boyer and Streng [@MR3376741], and Quinn and Verjovsky [@MR4091535].


Having access to the algorithm in this package, which is valid for any totally real number field, 
opens up for several new research directions and generalisations of previous research. 
Some of the direct applications to be pursued by the package author and collaborators 
lie in the field of explicit formulas and computational aspects of non-holomorphic Hilbert modular forms. 

# Implementation

The package `hilbertmodgroup` is mainly written in Python with some parts in Cython [@behnel2011cython]. 
It is intended to run as a package inside SageMath [@sage] as it 
makes heavy use SageMath's implementation of 
number fields, which is in turn is in many cases using the backend from PARI/gp [@PARI2].

# Documentation and Examples
All functions are documented using docstrings with integrated doctests 
following the guide for SageMath development. 
In addition, the directory `/examples` contains  
Jupyter notebooks illustrating the use of the package with 
a selection of fundamental examples, corresponding to examples presented in [@reductionalghilbert].

# Acknowledgements
This work was partially supported by the Engineering and Physical Sciences Research Council [EP/N007360/1] and [EP/V026321/1].

# References