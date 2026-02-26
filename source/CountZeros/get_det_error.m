function [f_error,df_error] = get_det_error(pulse4D,Eu)

    h1coeffs = [Eu.U_1_cheb(1,1);2*Eu.U_1_cheb(2:end,1)];
    h2coeffs = [Eu.U_1_cheb(1,2);2*Eu.U_1_cheb(2:end,2)];
    h4coeffs = [Eu.U_1_cheb(1,4);2*Eu.U_1_cheb(2:end,4)];
    phiP1coeffs = [Eu.U_vpp_cheb(1,1);2*Eu.U_vpp_cheb(2:end,1)];
    phiP2coeffs = [Eu.U_vpp_cheb(1,2);2*Eu.U_vpp_cheb(2:end,2)];
    phiP3coeffs = [Eu.U_vpp_cheb(1,3);2*Eu.U_vpp_cheb(2:end,3)];


    f_error = pulse4D.r*(norm(h1coeffs,1) + norm(h4coeffs,1)) ...
                    + Eu.r*(norm(phiP1coeffs,1) + norm(phiP2coeffs,1)) ...
                    + 2*pulse4D.r*Eu.r;
    df_error = pulse4D.r*(norm(h1coeffs,1) + norm(h2coeffs,1)) ...
                    + Eu.r*(norm(phiP1coeffs,1) + norm(phiP3coeffs,1)) ...
                    + 2*pulse4D.r*Eu.r;

end