function [output] = Cauchy_Mat_rep_trunc(Coeff_in,N_out)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% N - the dimension of the 

[n,m] = size(Coeff_in); % Assume square
N_inner=n;
N_inner=N_inner-1;
% Assume a and input is size N

toep_depth = N_out+1;
toep_width = N_out+1;


output = zeros((N_out+1)^2,(N_out+1)^2);

% Toeplitz_matrix_list = zeros()

% Make the a column of N, NxN matrixes given by Toepliz Matrices
toep_row_top = 0*(1:toep_depth);
toep_col = 0*(1:toep_depth);
for i = 1:N_inner+1
    toep_col(1:N_inner+1) = Coeff_in(1:N_inner+1,i);
    toep_row_top(1)=toep_col(1);
    Toeplitz_matrix = toeplitz(toep_col ,toep_row_top); % THis matrix isn't big enough


    output((i-1)*toep_depth+1:i*toep_depth,1:toep_width) = Toeplitz_matrix ;
end

% Copy this going down across rows, like a toepliz matrix 
Mat_column = output(1:toep_depth^2,1:toep_width);


for i = 2:N_out+1 
    N_shift_down = (i-1)*toep_depth;
    N_shift_right = (i-1)*toep_width;
    output(N_shift_down +1:end,N_shift_right+1:N_shift_right+toep_width) = Mat_column(1:end-N_shift_down,:);
end

 
end




function [coeff_out ] = changeRep_Vector2Spec( vec_in ,N,J)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

coeff_out=reshape(vec_in ,N,J);

end

function [vec_out] = changeRep_Spec2Vector( coeff_in  )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[N,J]=size(coeff_in );
N;
J;
vec_out=reshape(coeff_in,N*J,1);

end