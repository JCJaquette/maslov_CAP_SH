function [F_out] = F_op(a_vec,mu_s,vs1_shift,vs2_shift,g_1,params)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% a_in % assume a_{0,0} = 0; a_{2,0} = alpha_1;  a_{1,1} = alpha_2.

mu=params(1);

% Get the index;
[N,M]=size(vs1_shift);
if N~=M
    disp("WARNING NON SQUARE DIMENSION!")
end
N=N-1;



%  Make the input a grid.
a_vec=[0*a_vec(1);a_vec];
a_2index = reshape(a_vec,(N+1),(N+1));

a_2index(1,1)=0;

alpha_1=a_2index(3,1);
alpha_2=a_2index(2,2);

a_2index(3,1)=0;
a_2index(2,2)=0;

% Make stuff Larger
a_2index(end+1,:)=0*a_2index(end,:);
a_2index(:,end+1)=0*a_2index(:,end);

g_1(end+1,:)=0*g_1(end,:);
g_1(:,end+1)=0*g_1(:,end);

vs1_shift(end+1,:)=0*vs1_shift(end,:);
vs1_shift(:,end+1)=0*vs1_shift(:,end);

vs2_shift(end+1,:)=0*vs2_shift(end,:);
vs2_shift(:,end+1)=0*vs2_shift(:,end);
N=N+1;
% Evaluate

% Leading Operator
L_entry = 0*a_vec(1)*zeros(N+1,N+1);
for m=0:N
    for n=0:N
        xi = -(1+((m-1)*mu_s+n*mu_s')^2)^2-mu;
        if (m==2&&n==0)||(m==1&&n==1)
            xi=0*a_vec(1); % Going to be zero anywhay
        elseif m==1&&n==0
            xi = -1-mu+0*a_vec(1);
        end
        if isnan(xi)
            keyboard
        end
        L_entry(m+1,n+1)=xi;
    end
end


% g_1_w = conv2(g_1,a_2index);
g_1_w = Cauchyfft2(g_1,a_2index);
Initial_Sum = L_entry.*a_2index +vs1_shift*alpha_1 + vs2_shift*alpha_2+g_1; 

F_out = g_1_w;
F_out(1:N+1,1:N+1) = F_out(1:N+1,1:N+1)  +Initial_Sum;
 
end
