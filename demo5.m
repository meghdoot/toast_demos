% Demo 5: compute data and Jacobian for CW fDOT

clc,close all,clear all,
addpath(genpath('../v2/'))
addpath(genpath('./'))

ths = 0;thd = (pi/4:pi/4:pi+3*pi/4)';
sp = [25*cos(ths) 25*sin(ths)];ns = size(sp,1);
mp = [25*cos(thd) 25*sin(thd)];nd = size(mp,1);

% create mesh
[vtx,idx,eltp] = mkcircle(25,6,32,2); % radius, sectors, rings, bdd
% create mesh object
mesh = toastMesh (vtx,idx,eltp);
mesh.SetQM(sp,mp);
qvec = mesh.Qvec('Neumann','Gaussian',2);
mvec = mesh.Mvec('Gaussian',2,1.4);
n = mesh.NodeCount();
geom.n = n;

% create inclusions
mua = 0.01*ones(n,1);
mus = 1*ones(n,1);
h = 0*ones(n,1);
ref = 1.4*ones(n,1);

[vtx,idx,~] = mesh.Data;
r= 4; % Radius
cx = (-12.5)*cosd(45); % Center Xcoordi
cy = (-12.5)*sind(45); % Center ycoordi sind(0+ii*5)
% Inhomogeneity Formation
Index1=find(sqrt((cx-vtx(:,1)).^2+(cy-vtx(:,2)).^2)<r);
h(Index1) = 1;

figure,mesh.Display(h);hold on;
plot(sp(:,1),sp(:,2),'r*'),plot(mp(:,1),mp(:,2),'b*')
title('Target h')

S = dotSysmat(mesh, mua, mus, ref, 0);
% --- compute excitation data -----
phi_e = S \ qvec;
% --- compute emission data -----
Src = phi_e .* h;
phi_fl = S \ Src;
fdata = mvec.' * phi_fl;

figure,
plot(fdata);title('Measured fluorescence'),xlabel('Detector index')
