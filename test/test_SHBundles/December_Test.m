% K^3u*v+4 K^2u*Kv+Ku*(6K^2 v+2 v)+4u*(K^3v+K*v)

% v *K^3u + 4Kv *K^2u + (6K^2 v+2v) *Ku +4(K^3v+K*v) *u

% clear
N=15;
BOOL_PLOT = 0; 
delta_grid =0.1;

try
    lambda_s_1 =eigenvalues.s(:,1);
    lambda_s_2 =eigenvalues.s(:,2);
catch
    lambda_s_1 = -.1111 +1i*1.0062;
    lambda_s_2 =conj(lambda_s_1);
end


try
    N=max([N,length(v_coeff_pm)-1]);
end



tic


K = zeros(N+1,N+1);
weights = zeros(N+1,N+1);
row=(0:N)*lambda_s_1;
col=(0:N)'*lambda_s_2;
for i=1:N+1
    K(i,:)=K(i,:)+row;
    K(:,i)=K(:,i)+col;

    weights(i,:)=weights(i,:)+abs(row);
    weights(:,i)=weights(:,i)+abs(col);
end

% Define W_s
W_coeff = 1*(rand(N+1,N+1)-.5)+1i*(rand(N+1,N+1)-.5);

% weights=((1:N+1)'*(1:N+1));
weights(1,1)=1;
weights=5.^weights;

W_coeff =W_coeff ./weights;
W_coeff(1,1)=1;




try
    m=length(v_coeff_pm);
    W_coeff = zeros(N+1,N+1);
    W_coeff(1:m,1:m)=v_coeff_pm;
end


if BOOL_PLOT 
    figure
    val = plot2dTaylor( W_coeff  );
    title('W')
    % figure
    % x = -1:delta_grid :(1+1e-10);
    % y = -1:delta_grid :(1+1e-10);
    % [X,Y]=meshgrid(x,y);
    % surf(real(val)+log(abs(X)).*X.^3)
    

    figure
    surf(log(abs(W_coeff))/log(10));
    title('log coeff W')
    zlim([-16,1])
    colorbar
end



% return

% Make the inverse of W
% coeff_pad = zeros(4*(N+1),4*(N+1));
% coeff_pad(1:N+1,1:N+1)=W_coeff;

disp('Computing Inverse')
N_inv=1*(N+1);
W_coeff_inv = zeros(N_inv,N_inv); 
base =W_coeff(1,1);
try
    W_coeff_inv(1,1)=1/base;
catch
    W_coeff_inv(1,1)=1/mid(base);
end
% W_coeff_inv = ifft(fft(coeff_pad).^-1);
% W_coeff_inv = W_coeff_inv(1:N+1,1:N+1);





W_null = W_coeff;
W_null(1,1)=0;
tic

for alphaNorm=1:(N_inv-1)
    % alphaNorm
    % prod = Cauchyfft2(W_coeff_inv,W_null);
    prod = conv2(W_coeff_inv,W_null);
    for i=0:alphaNorm
        W_coeff_inv(1+i,alphaNorm-i+1) = -prod(1+i,alphaNorm-i+1)/base;
    end
end
toc

prod = conv2(W_coeff_inv,W_coeff);


% return
if BOOL_PLOT
    
    figure
    plot2dTaylor( W_coeff_inv  );
    title('W^{-1}')
    % figure
    % 
    % plot2dTaylor( W_coeff_inv, N);
    % figure
    % plot2dTaylor( conv2(W_coeff ,W_coeff_inv ), N);
    
    figure
    surf(log(abs(W_coeff_inv))/log(10));
    title('log coeff W^{-1}')
    zlim([-16,1])
    colorbar
end

% Initialize U
U_coeff = zeros(N+1,N+1);
U_coeff(1,1)=1; 



L_deriv = K.^3;
L_deriv_inv=L_deriv.^-1;
L_deriv_inv(1,1)=1;


W_coeff=intval(W_coeff);
U_coeff=intval(U_coeff);
K=intval(K);
toc

% initial output  vv = -( L e + Me ) = - Me 
vv = -M_operator( U_coeff, W_coeff,K,lambda_s_1);
vv=mid(vv);
Livv = Cauchyfft2(vv,W_coeff_inv);
Livv=mid(Livv);
Livv=Livv(1:N+1,1:N+1);

Livv =K.^(-3).*Livv;
Livv(1,1)=0;

Livv_vec = reshape(Livv,(N+1)*(N+1),1);


disp('Computing Bdagger')
Bdag = Bdagger( W_coeff,K,lambda_s_1,W_coeff_inv);

answer_vec  = mid(Bdag)\mid(Livv_vec(2:end));

answer_vec=[0;answer_vec  ];  

answer=reshape(answer_vec,N+1,N+1);
answer(1,1)=1;

% Compute Bounds
yy = (K.^3).*answer ;
yy(1,1)=0;
f_of_x_a = Cauchyfft2(W_coeff,yy);
f_of_x_b = + M_operator( answer, W_coeff,K,lambda_s_1);

f_of_x = f_of_x_a + f_of_x_b ;

sum(abs(f_of_x),'all')

return

tic
disp('Iterating')
for k=2:floor(N/2)
    disp(k)
    % disp('Continue')
    % pause
    % plot2dTaylor( W_coeff_square_trunc, N);
    M = -M_operator( U_coeff, W_coeff,K,lambda_s_1);
    

    next = Cauchyfft2(M,W_coeff_inv);
    % next = conv2(M,W_coeff_inv);
    next_trunc = next(1:N+1,1:N+1);
    next_trunc =L_deriv_inv.*next_trunc ;
    next_trunc(1,1)=1;
    for i=0:k
        U_coeff(1+i,k-i+1) = next_trunc(1+i,k-i+1);
    end
    % U_coeff=next_trunc;
    % U_coeff(1:8,1:8)
    
    out_final = Cauchyfft2(W_coeff, L_deriv.*U_coeff)+M_operator( U_coeff, W_coeff,K,lambda_s_1);
    % out_final = conv2(W_coeff, L_deriv.*U_coeff)+M_operator( U_coeff, W_coeff,K,lambda_s_1);
    result=sum(abs(out_final),'all');
    sup(result)

    % if k==3
    %     XXX=U_coeff;
    %     return
    % end
    
end
toc

N_large = 2*N+length(W_coeff_inv);
K_large = zeros(N_large,N_large);
row=(0:(N_large-1))*lambda_s_1;
col=(0:(N_large-1))'*lambda_s_2;
for i=1:N_large
    K_large(i,:)=K_large(i,:)+row;
    K_large(:,i)=K_large(:,i)+col;
end

L_deriv_large = K_large.^3+2*lambda_s_1*K_large.^2+2*(1+lambda_s_1^2)*K_large;
L_deriv_inv_large=L_deriv_large.^-1;
L_deriv_inv_large(1,1)=1;

% out_final = Cauchyfft2(W_coeff, L_deriv.*U_coeff)+M_operator( U_coeff, W_coeff,K,lambda_s_1);
out_final = conv2(W_coeff, L_deriv.*U_coeff)+M_operator( U_coeff, W_coeff,K,lambda_s_1);
sum(abs(out_final),'all')

% out_final_Y0=L_deriv_inv_large.*Cauchyfft2(out_final,W_coeff_inv);
out_final_Y0=L_deriv_inv_large.*conv2(out_final,W_coeff_inv);
sum(abs(out_final_Y0),'all')



if BOOL_PLOT
    figure
    plot2dTaylor( next_trunc  );
    title('New U')
    
    
    figure
    surf(log(abs(U_coeff))/log(10));
    colorbar
    title('Coeff of solution')


    
    figure
    surf(log(abs(out_final))/log(10));
    colorbar
    title('raw defect')



    figure
    surf(log(abs(out_final_Y0))/log(10));
    colorbar
    title('Y0 defect')

end

% W_coeff=U_coeff;
% 
% coeff_pad = zeros(4*(N+1),4*(N+1));
% coeff_pad(1:N+1,1:N+1)=W_coeff;
% 
% W_coeff_square = ifft(fft(coeff_pad).^2);
% W_coeff_square_trunc = W_coeff_square(1:2*N+1,1:2*N+1);
% figure
% plot2dTaylor( W_coeff_square_trunc, N);

toc


U_coeff_integrate=U_coeff./(K-2*lambda_s_1);
U_coeff_integrate(1,3)=0;

V_u0 = conv2(U_coeff_integrate,W_coeff); % mult by e^lambda_u
V_u0 = V_u0(1:N+1,1:N+1);
V_u1 = (K-lambda_s_1).*V_u0;
V_u2 = (K-lambda_s_1).*V_u1;
V_u3 = (K-lambda_s_1).*V_u2;


figure
val0 = plot2dTaylor( V_u0);
figure
val1 = plot2dTaylor( V_u1);
figure
val2 = plot2dTaylor( V_u2);
figure
val3 = plot2dTaylor( V_u3);

base=[V_u0(1,1);V_u1(1,1);V_u2(1,1);V_u3(1,1)];
base=base/norm(base);
base/base(1)
e_u=eigenvectors.u(:,1);
e_s=eigenvectors.s(:,1);
e_u=e_u/norm(e_u);
e_u/e_u(1)

normy = (abs(val0).^2+abs(val1).^2+abs(val2).^2+abs(val3).^2).^.5;
val0=val0./normy;
val1=val1./normy;
val2=val2./normy;
val3=val3./normy;

CombineVectors = zeros(length(val3),length(val3),4);
CombineVectors(:,:,1)=  val0;
CombineVectors(:,:,2)=  val1;
CombineVectors(:,:,3)=  val2;
CombineVectors(:,:,4)=  val3;

vecR = real(CombineVectors);
vecI = imag(CombineVectors);

vecR=vecR./(vecR(:,:,1).^2+vecR(:,:,2).^2+vecR(:,:,3).^2+vecR(:,:,4).^2).^.5;

dotty = 0*vecR(:,:,1);
for i=1:4
    
    dotty = dotty +vecR(:,:,i).*vecI(:,:,i);
    
end
sum(abs(dotty),'all')

for i=1:4
    vecI(:,:,i) = vecI(:,:,i) -dotty.*vecR(:,:,i);
end
vecI=vecI./(vecI(:,:,1).^2+vecI(:,:,2).^2+vecI(:,:,3).^2+vecI(:,:,4).^2).^.5;

dotty = 0*vecR(:,:,1);
for i=1:4
    dotty = dotty +vecR(:,:,i).*vecI(:,:,i);
end
sum(abs(dotty),'all')


averageR = reshape(mean(mean(vecR)),[4 1]);
averageI = reshape(mean(mean(vecI)),[4 1]);

varR=(vecR-mean(mean(vecR)));
varI=(vecI-mean(mean(vecI)));
V3=reshape(mean(mean( varI.*varR)),[4 1]);
V3=V3/norm(V3)

averageR=averageR/norm(averageR);

averageI = averageI-dot(averageI,averageR)*averageR;
averageI=averageI/norm(averageI);

% V3 = [0;0;1;0];
V3 = V3-dot(averageI,V3)*averageI;
V3 = V3-dot(averageR,V3)*averageR;
V3=V3/norm(V3)

V4 = [0;0;0;1];
V4 = V4-dot(averageI,V4)*averageI;
V4 = V4-dot(averageR,V4)*averageR;
V4 = V4-dot(V3,V4)*V3;
V4=V4/norm(V4);

AA = [ real(e_u)/norm(real(e_u)) imag(e_u)/norm(imag(e_u)) real(e_s)/norm(real(e_s)) imag(e_s)/norm(imag(e_s))];
% AA=AA/norm(AA);
AAi = inv(AA);
% AA= [averageR averageI V3 V4];
% transpose(AA)*AA;
vecR_new = 0*vecR;
vecI_new = 0*vecR;
mmm=length(vecR(:,1,1));
for i = 1:mmm
    for j = 1:mmm
        vecR_new(i,j,:)=AAi*reshape(vecR(i,j,:),[4,1]);
        vecI_new(i,j,:)=AAi*reshape(vecI(i,j,:),[4,1]);
    end
end
    

% delta_graph = .1;
x = -1:delta_grid:1;
y = -1:delta_grid :1;
[X,Y]=meshgrid(x,y);
figure
hold on
quiver3(X,Y,0*X,real(val0),real(val1),real(val3))
quiver3(X,Y,0*X,imag(val0),imag(val1),imag(val3))
hold off

% figure
% hold on
% quiver3(X,Y,0*X,vecR(:,:,1),vecR(:,:,2),vecR(:,:,4))
% quiver3(X,Y,0*X,vecI(:,:,1),vecI(:,:,2),vecI(:,:,4))
% hold off

figure
hold on
quiver3(X,Y,0*X,vecR_new(:,:,1),vecR_new(:,:,2),vecR_new(:,:,4))
quiver3(X,Y,0*X,vecI_new(:,:,1),vecI_new(:,:,2),vecI_new(:,:,4))
hold off

sum(abs(conv2(W_coeff_inv,K.*W_coeff)),'all')
sum(abs(K_large(1:2*N+1,1:2*N+1).*conv2(W_coeff_inv,K.*W_coeff)),'all')

return

function prod_trunc = Cauchyfft2( A,B)
N1 = length(A);
N2 = length(B);
N = 2*max(N1,N2);

coeff_pad = intval(0)*zeros(2*N,2*N);
A_pad = coeff_pad;
B_pad = coeff_pad;
A_pad(1:N1,1:N1)=A;
B_pad(1:N2,1:N2)=B;


% 
% prod = (fft(fft(A_pad,1).',1).') .* (fft(fft(B_pad,1).',1).');
% prod = fft(fft(prod,-1).',-1).';
prod = (verifyfft(verifyfft(A_pad,1).',1).') .* (verifyfft(verifyfft(B_pad,1).',1).');
prod = verifyfft(verifyfft(prod,-1).',-1).';
prod_trunc = prod(1:N,1:N);

end

function Bdag = Bdagger( W,K,lambda,wI)
N=length(W);
% 2N truncation
multiple = 1;
N_large = multiple*N;

Bdag = intval(0)*zeros(N_large^2-1,N_large ^2-1);

% h=zeros(N_large,N_large);

Bdag_data = intval(0)*zeros(N_large,N_large,N_large^2-1);

for i=1:N_large
    i
    for j = 1:N_large
        
        if i==1 && j==1 
            continue
        end

        h=zeros(N_large,N_large);
        h(i,j)=1;
    
        Mh = M_operator( h, W,K,lambda);
        wIMh = Cauchyfft2(wI,Mh);
        wIMh=wIMh(1:N,1:N);


        AMh = K.^(-3).*wIMh ;

        vec_out=reshape(AMh,N*N,1);

        Bdag_data(i,j,:)=vec_out(2:end);

        % index= (i-1)*N_large+(j-1);
        % Bdag(:,j)=vec_out(2:end);
    end
end

for i=1:N_large
    i
    for j = 1:N_large
        if i==1 && j==1 
            continue
        end
        index= (j-1)*N_large+(i-1);
        Bdag(:,index)=Bdag_data(i,j,:);
    end
end
Bdag =Bdag + eye(N_large^2-1);


end



function out = M_operator( U, W,K,lambda)
% U - w_u
% W - w_s
K_2 = K.^2;
K_3 = K.^3;
K_2(1,1)=intval(0);
K_3(1,1)=intval(0);
a2=(2*lambda*W+4*K.*W);
a1=((1+lambda^2)*W+ 4*lambda*K.*W+6*K_2.*W);
a0=4*( (1+lambda^2)*K.*W+K_3.*W);

B = Cauchyfft2(a2,K_2.*U);
C = Cauchyfft2(a1,K.*U);
D = Cauchyfft2(a0,U);

% U_coeff
% B = 4*conv2(K.*W,K.^2.*U);
% C = conv2(4*lambda*K.*W+6*K.^2.*W,K.*U);
% D = 4*conv2((1+lambda^2)*K.*W+K.^3.*W,U);

out = B+C+D;
end

function val = plot2dTaylor( coeff )
N=length(coeff)-1;
delta_grid = .1;
x = -1:delta_grid :1;
y = -1:delta_grid :1;
[X,Y]=meshgrid(x,y);
val = 0*Y'*Y;

X_powers = zeros(length(x),length(x),N+1);
Y_powers = zeros(length(x),length(x),N+1);
for i = 0:N
    X_powers(:,:,i+1)=X.^i;
    Y_powers(:,:,i+1)=Y.^i;
end

XtoN=1+0*X;

for i=0:N
    val_local = 0*val;
    for j=0:N
        val_local=val_local+coeff(i+1,j+1)*Y_powers(:,:,j+1);
    end
    val = val + val_local.*X_powers(:,:,i+1);
end

surf(X,Y,real(val))

end



