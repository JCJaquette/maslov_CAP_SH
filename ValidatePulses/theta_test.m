hold on
plot_manifold(mflds.unstable.coeffs, params.mfld.order, 'red')
plot_manifold(mflds.stable.coeffs, params.mfld.order, 'blue')
plot3(yo1 ,yo2 ,yo4 ,'LineWidth',1,'color','black')
plot3(right_endpt_s(1),right_endpt_s(2),right_endpt_s(4),'. black','MarkerSize',16);
plot3(left_endpt_u(1),left_endpt_u(2),left_endpt_u(4),'. black','MarkerSize',16);
theta_complex = params.rho*exp(1i*new_y.psi);
theta_phi1 = real(theta_complex);
theta_phi2 = imag(theta_complex);
ptpt = get_manifold_point(mflds.stable.coeffs,theta_phi1,theta_phi2,15);
ptpt = real(ptpt);
plot3(ptpt(1),ptpt(2),ptpt(4),'*')