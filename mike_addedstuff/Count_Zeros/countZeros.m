close all
clear
clc
global_min = -3;
global_max = 3; 
n = 500;
dom = global_min:1/n:global_max;
% f = @(x) x.^2 + sin(50*x)-2;
% df = @(x) 2*x + 50*cos(50*x);
e = 10.^([0:9]'); %#ok<NBRAK1>
% a = 10*rand(10,1)./e;
% A = IntegrateCheb(a);
A = [0;1.60489867823427;0.0366042541650144;0.0146799834264349;0.000936108606866844;8.31884488841465e-05;9.38366320848015e-07;5.67664670792784e-07;1.35106762535338e-08;3.50996476709789e-09;2.78005655285905e-10];
a = [3.29871715200074;0.153917367641776;0.0889197955322078;0.00750035098171809;0.000839894973598372;1.14821267833384e-05;8.01048475690674e-06;2.21730933162258e-07;6.31793658077621e-08;5.56011310571810e-09];
df = @(x) chebSum(a,x);
f = @(x) chebSum(A,x);

[roots,indeces] = FindZero1D(dom,f(dom));
grid_breakup = [global_min;roots;global_max];
rs = 0;
unique0 = [rs,rs];

for i = 1:length(roots)

    rs = RadiiPoly1D(f,df,intval(1)*roots(i),[grid_breakup(i), grid_breakup(i+2)]);
    unique0(i,:) = [roots(i) - rs.sup, roots(i) + rs.sup];
    plot(unique0(i,:),[0,0],'color','blue')

end

if 0*unique0 == [0,0]
    if unique0(1) <= dom(1)
        newDom = infsup(unique0(2),dom(2));
    elseif unique0(2) >= dom(2)
        newDom = infsup(dom(1),unique0(1));
    end
else
    newDom = infsup([global_min;unique0(:,2)],[unique0(:,1);global_max]);
end

stop = 0; rejects = [];

while isempty(newDom) == 0

    b = newDom(end);
    if abs(f(b)) > 0
        newDom(end) = [];
        stop = 0;
    else 
        newDom(end) = [];
        newDom = [newDom;
                  infsup(b.inf, (b.inf + b.sup)/2);
                  infsup((b.inf + b.sup)/2, b.sup)]; %#ok<AGROW>

        stop = stop + 1;
    end

    if stop == 20
        
        rejects = newDom(end);
        newDom(end) = [];

    end

if length(newDom) > 10000
    disp('something is wrong')
    return
end

end