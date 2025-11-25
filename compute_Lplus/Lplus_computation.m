clear

all_bundles %This gets the bundles and manifolds
clearvars -except mflds bndl

getICs_runNewton_validate %This gets the pulse
clearvars -except mflds bndl new_y

chebInt %This gets U_{\varphi'} and U_1
clearvars -except mflds bndl new_y phiPrime_cheb h_cheb 

% load('saved_stuff/sig0_1.mat'); %Pick sig0 depending on the n
% load('saved_stuff/sig0_2.mat');
load('saved_stuff/sig0_3.mat');

computeLplus(params,bndl,mflds,U_1,sig0);
