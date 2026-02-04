function psoln = get4Dpsoln(timevec,pulse)
% Get approximate derivatives with finite differences

psoln = [timevec,pulse];
dx = timevec(2) - timevec(1);
pulsei = pulse;
for i = 3:5
    pad = [0;0;pulsei;0;0];
    pulsei = (-pad(5:end) + 8*pad(4:end-1) - 8*pad(2:end-3) + pad(1:end-4))/(12*dx);
    psoln(:,i) = pulsei;
end

deriv2 = (pulse(3:end)-pulse(1:end-2))/(2*dx);

% figure 
% plot(timevec(2:end-1), deriv2 )
% hold on 
% plot(timevec(2:end-1), psoln(2:end-1,3) )

end