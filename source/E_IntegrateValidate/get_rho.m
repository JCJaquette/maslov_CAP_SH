function rho = get_rho(params,mflds,psoln,phi_cheb,nonzero)
%gets the error rho in the parameter lambda of corollary 3.1

maxphi = max(abs(psoln.phi1),abs(psoln.phi2)); 
mani_tail_error = mflds.unstable.r_min; 
mani_coeff_error = reshape(sum(sum(rad(mflds.unstable.coeffs))),[4,1]);
mani_error = max(mani_coeff_error + mani_tail_error*ones(4,1));%check this?

IC_error = 2*pi/log(1/maxphi) * mani_error;%Error of manifold's tangent bundle, lemma 8 in hexagon paper

psoln.r = 1.8e-11;
tail = psoln.a1' - [phi_cheb;zeros(params.Eu.order-nonzero,1)];
tail_error = vectorDelta1norm(tail,params.del);

pulse_error = psoln.r + tail_error;

rho = max(IC_error,pulse_error);

end