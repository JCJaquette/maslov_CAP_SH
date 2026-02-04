This directory contains the code that generates the main results, aiming to provide a computer assisted proof of stability for pulses to the Swift-Hohenberg via conjugate points. They are as follows:

* [all_bundles.m](all_bundles.m): Computes the (un)stable manifold, L_minus, and the stable/unstable bundles over the stable manifold.  Primarily relies on the code in [source/InvariantManifolds](source/InvariantManifolds) and [source/SwiftHohenbergBundles](source/SwiftHohenbergBundles). 

* [conjugate_points.m](conjugate_points.m): Computes (via standard numerics) the number of conjugate points and approximate unstable eigenvalues associated to three pulse solutions of the Swift-Hohenberg equation. Two of the solutions lie in the non-snaking parameter region and are unstable; one lies in the snaking parameter region and is spectrally stable. This script relies on the code contained in [Fourier Series](source/@FourierSeries), [Pulse Solution](source/@PulseSolution), and [Conjugate Points](source/@ConjugatePoints).

* [manifold_validation.m](manifold_validation.m): Computes a Taylor approximation for the stable and unstable manifolds of the Swift-Hohenberg equation using the parameterization method and bounds the error on these approximations via a computer assisted proof. Primarily relies on the code in [source/InvariantManifolds](source/InvariantManifolds).

* [pulse_validation.m](pulse_validation.m): Computes an approximation of a homoclinic orbit and bounds the error as discussed in Chapter 9. This relies on the code in [source/InvariantManifolds](source/InvariantManifolds) and [source/PulseValidation](source/PulseValidation).

* nu_1p6_mu_0p2 - Contains the output of all_bundles, including timing information. 
