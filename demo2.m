% Demo 2: data generation, in-homogeneous properties

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
ref= 1.4*ones(n,1);

% Specify inhomogenities/inclusions in domain
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

% Specify optical properties
mua0 = 0.01;mus0 = 1;kap0 = 1/(3*(mua0+mus0));
mua = mua0*ones(n,1);kap = kap0*ones(n,1);
mua(Index1) = mua0*20;kap(Index2) = kap0*20; % specifying inclusions
mus = (1/3)*(1./kap-3*mua);

% Solve forward problem
K = dotSysmat(mesh,mua,mus,ref,100);
Phi1 = K\qvec; % Fluence
Gamma1 = mvec.' * Phi1; % Exitance
data1 = [log(abs(Gamma1(:)));angle(Gamma1(:))];

% add noise
figure, subplot(121),mesh.Display(mua);title('target \mu_a');
subplot(122),mesh.Display(kap);colorbar,title('target \kappa');

figure,
subplot(121),mesh.Display(log(abs(Phi1(:,1)))),title('LAmp')% log-amplitude
subplot(122),mesh.Display(angle(Phi1(:,1))),title('Phs')% phase

figure,
plot([log(abs(Gamma1(:)));angle(Gamma1(:))]),title('Data')

