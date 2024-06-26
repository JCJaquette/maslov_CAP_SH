This directory contains the code that generates the main results discussed in this dissertation. They are as follow: 

* [bistable_vector_bundles.mlx](bistable_vector_bundles.mlx): Computes the invariant vector bundles for the bistable equation using the reduction of order method discussed in Chapter 11.2 and compares these solutions to the known analytic expressions. This script primarily relies on the code contained in [source/BistableBundles](source/BistableBundles).

* [conjugate_points.m](conjugate_points.m): Computes the number of conjugate points and approximate unstable eigenvalues associated to three pulse solutions of the Swift-Hohenberg equation. Two of the solutions lie in the non-snaking parameter region and are unstable; one lies in the snaking parameter region and is spectrally stable. This script relies on the code contained in [Fourier Series](source/@FourierSeries), [Pulse Solution](source/@PulseSolution), and [Conjugate Points](source/@ConjugatePoints).

* [L_minus.m](L_minus.m): Computes an estimate on $L_-$ for the three pulse solutions, where $L_-$ is the value for which there are no conjugate points for $x < -L_-$. This script relies on the code contained in [source/L_minus](source/L_minus).

* [manifold_validation.m](manifold_validation.m): Computes a Taylor approximation for the stable and unstable manifolds of the Swift-Hohenberg equation using the parameterization method and bounds the error on these approximations via a computer assisted proof. This is as discussed in Chapter 9.1 and primarily relies on the code in [source/InvariantManifolds](source/InvariantManifolds).

* [pulse_validation.m](pulse_validation.m): Computes an approximation of a homoclinic orbit and bounds the error as discussed in Chapter 9. This relies on the code in [source/InvariantManifolds](source/InvariantManifolds) and [source/PulseValidation](source/PulseValidation).

* [swift_hohenberg_vector_bundles.m](swift_hohenberg_vector_bundles.m): Computes the non-resonant vector bundles for the Swift-Hohenberg equation and implements the necessary functions for computing the resonant vector bundles using the reduction of order method. This is as disucssed in Chapter 11.3 and relies on the code in [source/SwiftHohenbergBundles](source/SwiftHohenbergBundles). 


