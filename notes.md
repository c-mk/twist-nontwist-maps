# twist-nontwist-maps


# task

reproduce the results on the twist and non-twist maps


# PAPER NOTES

**title: PJ Morrison 2000 - Magnetic field lines, Hamiltonian dynamics, and nontwist systems**

1. Magnetic field lines do not simply close or go to infinity. 
2. Magnetic field lines are in fact a Hamiltonian system with chaotic trajectories. 
   3. because the map from the p - thetha plane to itself is an area preserving map which is identical to the Hamiltonian trajectory with the role of time being played by the azithmul coordinate phi
3. Hamiltonian systems with two degrees of freedom have a position q and a momenta p have 2 coordinates (p1, q1), (p2, q2)
4. Taking a 2d slice of this trajectory is a Poincare section
5. Because irrational numbers infinitely outnumber rational numbers in mathematics, a randomly selected or "typical" field line will have an irrational rotational transform. Instead of closing, it will endlessly trace out a dense web that eventually covers an entire 2D surface (known as a KAM torus) inside the volume

**useful definitions**

- zero dimensional space - a single point
- one dimensional space - a straight line 
- two dimensional space - covers a plan in two directions (x, y)
- three dimensional space - covers length, depth (front to back length) vs width (left to right width), height (top to bottom)
- four dimensional space (x,y,z,t) time!
- euclidean n-dimensional space - any number of dimensions 
- vector spaces - finite number of dimensions
- hilbert spaces - infinite dimensional spaces used in quantum mechanics
- manifolds - spaces that look like euclidean space near every point but have a different shape
- phase space describes the complete set of variables nedded to predict the future behavior of a dynamic system
- closed magnetic field line - rational ratio (toroidal to poloidal (long, short)) of rotational transform indicative of the field line closing in on itself
- open magnetic field line - irrational ratio of rotational transform indicative of the magnetic field line spiriling forever and never perfectly overlapping itself
- divergence-free magnetic field line - a field like that has no beginning or end (specifically Gauss's Law for magnetism) 


**summary**
In a toroidal magnetized plasma (tokamak, stellarator), magnetic field lines can be described as trajectories of a Hamiltonian system 

— following a field line around the torus is equivalent to iterating an area-preserving map (a Poincaré map) in a 2D cross-section. This comes directly from the fact that field-line flow conserves magnetic flux, so it's symplectic/area-preserving, just like a Hamiltonian phase-space flow.

Greene's Criterion (Greene, 1979) is a method from KAM (Kolmogorov-Arnold-Moser) theory used to determine whether an invariant torus (a smooth, continuous magnetic flux surface) survives or breaks apart into a fractal remnant called a cantorus.

**questions:**

what is the jacobian and its relation with plasma physics?
is the residue of Greene's Criterion a magnetized plasma left : it is not a physical goo, instead it is a measurement of stability of a specific periodic orbit (elliptic vs hyperbolic)
what is meant by the phase space action principle that derives the Hamiltonian's equations?

The term “non-resonant divertor” was coined in 2015 by Boozer [3]. The definition emphasised (1) the relationship between sharp curvature of flux surfaces and X-points in the magnetic field and (2) the importance of making perturbations to the magnetic field which are not resonant with the rotational transform of the field, and therefore do not produce an island divertor