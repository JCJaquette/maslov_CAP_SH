function [roots,flag] = FindZero1D(dom,f)
%take in f values, check if there is zero in domain

plot(dom,f,'Color','black')
hold on

Zs = []; 

for i = 1:length(f)-1 %does f change sign

    if f(i)*f(i+1)<0 || abs(f(i)) < 10^-8
        Zs = [Zs;[dom(i),i]];
        plot(Zs(end,1),0,'.','Color',"green")
    end

end

if isempty(Zs) == 0
    roots = Zs(:,1);
    flag = 0;
else
    roots = 0;
    flag = 1;
    disp('no roots found')
end



end