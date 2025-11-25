
params = setparams();

mflds = get_mflds(params);

hold on

[u_pts,u_p1p2s] = plot_manifold(mflds.stable.coeffs,25,'blue');
s_pts = plot_manifold(mflds.unstable.coeffs,25,'red');

%plot3(pulseSoln.normalForm.sol(:,1),pulseSoln.normalForm.sol(:,2),pulseSoln.normalForm.sol(:,4),'black');
%plot3(pulseSoln.normalForm.sol(:,1),-pulseSoln.normalForm.sol(:,2),-pulseSoln.normalForm.sol(:,4),'black');


%%

normsu = zeros(1,59);
for i = 1:59
c1 = reshape(u_pts(60,i,:),[4,1]);
c2 = reshape(u_pts(60,i+1,:),[4,1]);
normsu = norm(c1-c2);
end

boundary_distance= mean(normsu);
boundary_distance=boundary_distance*.6;

n = length(pulseSoln.normalForm.sol(:,1));
d_u = zeros(n,1);
pts_u = zeros(n,4);
d_s = zeros(n,1);
pts_s = zeros(n,4);
mps = zeros(n,2);

for k = 1:n  
    [pts_u(k,:), d_u(k), mps(k,:)] = closestpt(u_pts,(pulseSoln.normalForm.sol(k,:))');
end
for k = 1:n
    if d_u(k) < boundary_distance
        break
    end
end

L_pl = pulseSoln.normalForm.time(k);
pt_u = pts_u(k,:);
Lsoln = pulseSoln.normalForm.sol(1:k,:);
phi1 = u_p1p2s(mps(k,1),mps(k,2),1);
phi2 = u_p1p2s(mps(k,1),mps(k,2),2);
psi = k/60 * 2*pi;


pulseSoln.normalForm.sol(:,2) = -pulseSoln.normalForm.sol(:,2);
pulseSoln.normalForm.sol(:,4) = -pulseSoln.normalForm.sol(:,4);

for k = 1:n  
    [pts_s(k,:), d_s(k)] = closestpt(s_pts,(pulseSoln.normalForm.sol(k,:))');
end
for k = 1:n
    if d_s(k) < boundary_distance
        break
    end
end

L_mn = -pulseSoln.normalForm.time(k);
pt_s = pts_s(k,:);
Lsoln = [pulseSoln.normalForm.sol(k:-1:1,:);
                    Lsoln];

pulseSoln.normalForm.sol(:,2) = -pulseSoln.normalForm.sol(:,2);
pulseSoln.normalForm.sol(:,4) = -pulseSoln.normalForm.sol(:,4);

hold on
u_pts = plot_manifold(mflds.stable.coeffs,25,'blue');
s_pts = plot_manifold(mflds.unstable.coeffs,25,'red');
plot3(Lsoln(:,1),Lsoln(:,2),Lsoln(:,4))

plot3(pt_s(1),pt_s(2),pt_s(4),'* black','MarkerSize',20);
plot3(pt_u(1),pt_u(2),pt_u(4),'* black','MarkerSize',20);

plot3(Lsoln(1,1),Lsoln(1,2),Lsoln(1,4),'* black','MarkerSize',20)
plot3(Lsoln(end,1),Lsoln(end,2),Lsoln(end,4),'* black','MarkerSize',20)





