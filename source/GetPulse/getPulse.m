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
S.fourier.M = 1000; 
S.fourier.tol = 1e-14; 
S.fourier.order = 500; 
S.time = 150; 

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

% perform Newton's method
S = BKNormalForm4d_halfline(S);


fulltime = [-flip(S.normalForm.time); S.normalForm.time(2:end)];
psoln = [fulltime,Newton_halfline(S)];
psoln = get4Dpsoln(psoln(:,1),psoln(:,2));



end