%pulse

pls1 = chebcoeff_to_function(new_x.a1);
pls2 = chebcoeff_to_function(new_x.a2);
pls3 = chebcoeff_to_function(new_x.a3);
pls4 = chebcoeff_to_function(new_x.a4);

hold on
plot3(pls1,pls2,pls4,'black');
hold on
plot3(pls1(1),pls2(1),pls4(1),"* red");
plot3(pls1(end),pls2(end),pls4(end),"* blue");

%% manifolds

u_maniTCs = calc_proj_coeff(evals.u,evecs.u,params);
plot_manifold(u_maniTCs,25,'red');

s_maniTCs = calc_proj_coeff(evals.s,evecs.s,params);
plot_manifold(s_maniTCs,25,'blue');





