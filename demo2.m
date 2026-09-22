clc, clear all, close all,

% create mesh
[vtx,idx,eltp] = mkcircle(25,6,32,2); % radius, sectors, rings, bdd
% create mesh object
mesh = toastMesh (vtx,idx,eltp);
% specify sources and detectors
ths = (-pi:2*pi/8:pi)';thd = 0.4-ths; 
sp = [25*cos(ths) 25*sin(ths)];mp = [25*cos(thd) 25*sin(thd)];

% FEM
mesh.SetQM(sp,mp);
qvec = mesh.Qvec('Neumann','Gaussian',2);
mvec = mesh.Mvec('Gaussian',2,1.4);    
n = mesh.NodeCount();
mua = 0.01*ones(n,1);    
mus= 1*ones(n,1);
ref= 1.4*ones(n,1);
K = dotSysmat(mesh,mua,mus,ref,100);

% solve Phi
Phi0 = K\qvec; % Fluence
Gamma0 = mvec.' * Phi0; % Exitance
data0 = [log(abs(Gamma0(:)));angle(Gamma0(:))];

% specify inhomogenities/inclusions in domain
[vtx,~,~] = mesh.Data;
% Circular Inhomogeneity Size
r= 4; % Radius
cx= 12.5; % Center Xcoordi
cy= 0; % Center ycoordi
% Inhomogeneity Formation
Index1=find(sqrt((cx-vtx(:,1)).^2+(cy-vtx(:,2)).^2)<r);
% Circular Inhomogeneity Size
r= 4; % Radius
cx= 0; % Center Xcoordi
cy= 12.5; % Center ycoordi
% Inhomogeneity Formation
Index2=find(sqrt((cx-vtx(:,1)).^2+(cy-vtx(:,2)).^2)<r);

mua0 = 0.01;mus0 = 1;kap0 = 1/(3*(mua0+mus0));
mua = mua0*ones(n,1);kap = kap0*ones(n,1);
mua(Index1) = mua0*2;kap(Index2) = kap0*2; % specifying inclusions
mus = (1/3)*(1./kap-3*mua);

K = dotSysmat(mesh,mua,mus,ref,100);
% solve Phi
Phi1 = K\qvec; % Fluence
Gamma1 = mvec.' * Phi1; % Exitance
data1 = [log(abs(Gamma1(:)));angle(Gamma1(:))];

deltay = data1-data0;

% add noise
deltay = deltay + 0.01*abs(deltay).*randn(length(deltay),1);

figure,
plot(deltay)

figure, subplot(121),mesh.Display(mua-mua0);title('target \delta\mu_a');
subplot(122),mesh.Display(kap-kap0);colorbar,title('target \delta\kappa');

%% Inverse problem

stdn = 0.01*abs(deltay); % assumed standard deviation of noise

Gamma_e = diag((stdn).^2);
iGamma_e = inv(Gamma_e);

mua = mua0*ones(n,1);mus = mus0*ones(n,1);
J = toastJacobian (mesh, [], qvec, mvec, mua, mus, ref, 100);

% White noise prior
Gamma_mua = ((0.02-0.01)/3)^2*eye(n);
Gamma_kap = ((0.66-0.33)/3)^2*eye(n);
Gamma_x = [Gamma_mua sparse(n,n);sparse(n,n) Gamma_kap];
iGamma_x = inv(Gamma_x);

x = (J'*iGamma_e*J + iGamma_x)\(J'*iGamma_e*deltay);

% reconstructions
muarecon = x(1:n);kaprecon = x(n+1:end);
figure, subplot(121),mesh.Display(muarecon);title('reco \delta\mu_a WN');
subplot(122),mesh.Display(kaprecon);colorbar,title('reco \delta\kappa WN');

[vtx,elem,~] = mesh.Data;
meshD = struct;
meshD.g = vtx;
r1 = 8; % correlation length set as approximately the size of target
Gamma_mua = PriorOU(meshD,(0.02-0.01)/3,r1);
Gamma_kap = PriorOU(meshD,(0.66-0.33)/3,r1);
Gamma_x = [Gamma_mua sparse(n,n);sparse(n,n) Gamma_kap];
iGamma_x = inv(Gamma_x);

x = (J'*iGamma_e*J + iGamma_x)\(J'*iGamma_e*deltay);

% reconstructions
muarecon = x(1:n);kaprecon = x(n+1:end);
figure, subplot(121),mesh.Display(muarecon);title('reco \delta\mu_a OU');
subplot(122),mesh.Display(kaprecon);colorbar,title('reco \delta\kappa OU');


