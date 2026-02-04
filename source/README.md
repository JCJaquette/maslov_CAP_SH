This directory contains the source code that generates the results presented [here](results). The
first 3 modules are as follows:  
* [Fourier Series](source/@FourierSeries): contains code for a `FourierSeries` object 
and operations.

* [Pulse Solution](source/@PulseSolution): contains code for a `PulseSolution`
object. This approximates the pulse solution to the Swift-Hohenberg equation 
using a Fourier series and Newton's method.

* [Conjugate Points](source/@ConjugatePoints): contains code for a `ConjugatePoints` 
object. This computes the conjugate points associated to the pulse solution. 


The following modules contain code that accompanies Part 2 of the dissertation. 
  
* [Invariant Manifolds](source/InvariantManifolds): contains code to compute the invariant stable and unstable manifolds for the Swift-Hohenberg equation.
  
* [L minus](source/Lminus): contains code to compute the bound on $L_-$.
  
* [Pulse Validation](source/PulseValidation): contains code to compute a stationary pulse solution of the Swift-Hohenberg equation and a computer assisted proof of its existence.
  
* [Sequence Space](source/SequenceSpace): contains helper functions for operations in sequence space.
  
* [Swift-Hohenberg Bundles](source/SwiftHohenbergBundles): contains code to compute the resonant and non-resonant solutions of the Swift-Hohenberg equation. 
  
* [Vector Field](source/VectorField): contains helper functions for expressing the Swift-Hohenberg equation, its vector field, and linearization about the origin. 
