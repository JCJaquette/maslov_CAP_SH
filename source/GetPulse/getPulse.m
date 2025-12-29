function psoln = getPulse(params)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% normal form branch -- set to 0, pi
S.normalForm.branch = params.xi; 

% vector field parameters
S.vfParams.nu = params.nu;
S.vfParams.mu = params.mu;
S.vfParams.lambda = 0; 

% fourier approximation parameters
S.fourier.M = 1500; 
S.fourier.tol = 1e-14; 
S.fourier.order = 500; 
S.time = 100; 

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(['We consider the pulse for parameter values nu=', ...
    num2str(params.nu), ', mu=', ...
    num2str(params.mu), ', and branch phi=', ...
    num2str(params.xi), '.'])
disp(' ')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                PULSE SOLUTION APPROXIMATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

disp('Performing Newtons method to obtain a Fourier approximation of the pulse solution.')
disp(' ')

% Get seed solution from the Burke Knobloch normal form
S = BKNormalForm4d_halfline(S);
% perform Newton's method

fulltime = [-flip(S.normalForm.time); S.normalForm.time(2:end)];
S= Newton_halfline(S);

full_sol = getDFunctionFromFourierCoeffs(S,S.fourier.full_coeff_from_half_newton, "full");

psoln = [fulltime,S.full_uout ];
psoln = get4Dpsoln(psoln(:,1),psoln(:,2));



end