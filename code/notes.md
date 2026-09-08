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

The term “non-resonant divertor” was coined in 2015 by Boozer. The definition emphasised (1) the relationship between sharp curvature of flux surfaces and X-points in the magnetic field and (2) the importance of making perturbations to the magnetic field which are not resonant with the rotational transform of the field, and therefore do not produce an island divertor


**questions:**

what is the jacobian and its relation with plasma physics?
is the residue of Greene's Criterion a magnetized plasma left : it is not a physical goo, instead it is a measurement of stability of a specific periodic orbit (elliptic vs hyperbolic)

what is meant by the phase space action principle that derives the Hamiltonian's equations?


# Supplemental Chapter in Plasma Confinement on Nonlinearity Processes


1. Magnetic field lines can be treated as a Hamiltonian system. When there are nested toroidal flux surfaces the field line trajectories are intrgrable, and the Hamiltonian corresponds to the toroidal flux. 

Most pertubabtion of an integral system will destroy integrability. The primary effect is the occurence of magnetic islands at every rational surface. When island sizes in width become to large compared to their speratrions, chaos arises.

2. When a Hamiltonian system is far from the integrable limit, statisical approaches can be used to treat chaotic trajectories. The simplest approach is analogous to quasilinear theory and leds to a diffusion coefficient depending only on the mean square level of the pertubing forces.

KAM theorem states that if X is smooth enough and small enough most of the original magnetic surfaces surve the pertubation. these surfaces are far from any of the rational surfaces. Every rational magnetic surface is destroted and replaced by an island. 

As pertubation approaches 0 most pertubations are perserved.

# more questions

What is island overlap? When the two islands are sperated by a distance much larger than the sum of their half widths?

Why do nonlinear effects distort island shapes and give rise to islands with other helicities (than single) ?

Are we thinking of the plasma as fluid for MHD?

What is island helicity?

(When there are multiple helicities, integration breaks down, why?)
# aug 29

Starlite: a stellerator designed to experiementally validate non-resonant divertors
https://arxiv.org/pdf/2603.18265

Divertor topology and vacuum vessel design for stellerators
https://arxiv.org/pdf/2607.27127


# sept 7

review papers

X-points are hyperbolic fixed points of the field-
line map; their stable and unstable invariant manifolds form the separatrix and the divertor legs that carry heat to the wall, and are
computed with SIMSOPT and pyoculus. 

Greene’s residue (spanning −0.19 to −0.35 across the three configurations) predicts how
quickly field lines travel along those manifolds, and in the high-ι case the manifolds fail to close, producing a homoclinic tangle and
edge chaos. 

Heat deposition is modeled with EMC3-Lite anisotropic transport and a Bohm sheath boundary condition, with 30 kW of exhausted power matching the two 15 kW, 2.45 GHz magnetrons. 

Student team tasks include reproducing the design A analysis from
the published SIMSOPT scripts, extending sensitivity studies to as-built coil geometry recovered from Project 1, and 
verifying strike
patterns and connection lengths for the selected vessel.

Skills to learn: Python scientific computing, dynamical systems and chaos, magnetic confinement theory, stellarator optimization, high-
performance computing.

**Timeline and Deliverables**

Project III: verify X-point locations from field measurement, 
divertor footprints measured and compared with
EMC3-Lite predictions, 
study effect of field error of the X-point locations. 


Each team will submit a 3-page written technical report,
calibrated hardware or documented code, 
and a presentation/poster suitable for an APS-DPP 
or IEEE PPPS undergraduate session.

# sept 8
submitted interest