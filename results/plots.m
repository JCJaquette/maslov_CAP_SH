figure
hold on
plot_manifold(mid(mflds.stable.coeffs),40,'blue',1)
plot_manifold(mid(mflds.unstable.coeffs),40,'red',1)
t = linspace(-1,1,200);  
for i = 1:length(t)
    x(i) = chebSum(pulse4D.a1, t(i));
    y(i) = chebSum(pulse4D.a2, t(i));
    z(i) = chebSum(pulse4D.a4, t(i));
end
plot3(x, y, z, 'black','LineWidth',1)
camorbit(-10,-15)
saveFigure

%%

t = linspace(-1,1,200);  
for i = 1:length(t)
    x(i) = chebSum(pulse4D.a1, t(i));
end

plot(t,x)
saveFigure

%%

t = linspace(-1,1,200); 
for i = 1:4
    for j = 1:length(t)
        U_vpp(j,i) = chebSum(Eu.U_vpp_cheb(:,i),t(j));
    end
end

for i = 1:4
    for j = 1:length(t)
        U_1(j,i) = chebSum(Eu.U_1_cheb(:,i),t(j));
    end
end

figure
hold on
tiledlayout(4,2)  

n = size(U_1,1);
x = linspace(-1,1,n);  

for i = 1:4

    nexttile
    plot(x, U_1(:,i), 'LineWidth', 1)
    title(['$U\_1$ component ', num2str(i)],'Interpreter', 'latex')
    grid on

    nexttile
    plot(x, U_vpp(:,i), 'LineWidth', 1)
    title(['$U\_{\varphi\prime}$ component ', num2str(i)],'Interpreter', 'latex')
    grid on
end

saveFigure