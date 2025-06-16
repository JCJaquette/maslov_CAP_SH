close all
clear

global_min = -1;
global_max = 1; 
n = 500; 
dom = global_min:1/n:global_max;
domINT = infsup(global_min,global_max);

% f = @(x) x.^2 + sin(20*x);
% df = @(x) 2*x + 20*cos(20*x);
% f_error = 0; df_error = 0;

% seriesLength = 30;
% sgns = (-1).^round(rand(seriesLength,1)); 
% a = 10*rand(seriesLength,1);
% a = a.*sgns;
% a(1) = 0;
% A = IntegrateCheb(a);
% 
% df = @(x) chebSum(a,x);
% f = @(x) chebSum(A,x);

[params,h_cheb,phi_cheb,phiPrime_cheb] = getparamsCount(1);

f_error = infsup(-params.f_error,params.f_error);
df_error = infsup(-params.df_error,params.df_error);

A_cheb = chebstar2fft(h_cheb(:,1),phiPrime_cheb(:,4)) ...
    - chebstar2fft(h_cheb(:,4),phiPrime_cheb(:,1));
APrime_cheb = chebstar2fft(h_cheb(:,1),phiPrime_cheb(:,2)) ...
    - chebstar2fft(h_cheb(:,2),phiPrime_cheb(:,1));

f = @(x) chebSum(A_cheb,x);
df = @(x) chebSum(APrime_cheb,x);

[roots,flag] = FindZero1D(dom,f(dom));
grid_breakup = [global_min;roots;global_max];
rs = 0;
unique0 = [rs,rs];

for i = 1:length(roots)

    rs = RadiiPoly1D(@(x) f(x) + f_error,@(x) df(x) + df_error, ...
        intval(1)*roots(i),[grid_breakup(i), grid_breakup(i+2)]);
    unique0(i,:) = [roots(i) - rs.sup, roots(i) + rs.sup];
    plot(unique0(i,:),[0,0],'color','blue')

end

if flag == 0
    zeroCount = length(unique0(:,1));
else
    zeroCount = 0;
end
unique0 = infsup(unique0(:,1),unique0(:,2));

newDom = INTsetminus(domINT,unique0);

rejects = intval(0)*[];

tol = 1e-5;
while isempty(newDom) == 0

    b = newDom(end);
    if abs(f(b)) > 0
        newDom(end) = [];
    else 
        newDom(end) = [];
        newDom = [newDom;
                  infsup(b.inf, (b.inf + b.sup)/2);
                  infsup((b.inf + b.sup)/2, b.sup)];                        %#ok<AGROW>
    end
    if rad(b)<tol
        rejects(end+1) = newDom(end);
        newDom(end) = [];
    end

end

if isempty(rejects) == 1
    disp('all roots on domain found')
    zeroCount
end