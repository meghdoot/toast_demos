% Demo 3: time-domain data generation

clc, clear all, close all,

% create mesh
[vtx,idx,eltp] = mkcircle(25,6,32,2); % radius, sectors, rings, bdd
% create mesh object
mesh = toastMesh (vtx,idx,eltp);
% specify sources and detectors
ths = 0;thd = pi; 
sp = [25*cos(ths) 25*sin(ths)];mp = [25*cos(thd) 25*sin(thd)];

% FEM
mesh.SetQM(sp,mp);
qvec = mesh.Qvec('Neumann','Gaussian',2);
mvec = mesh.Mvec('Gaussian',2,1.4);    
n = mesh.NodeCount();
ref= 1.4*ones(n,1);

% Specify optical properties
mua0 = 0.01;mus0 = 1;
mua = mua0*ones(n,1);mus = mus0*ones(n,1);

% Obtain FE matrices
A = dotSysmat(mesh,mua,mus,ref,0);
M = mesh.Massmat;

dt = 10;      % step size [ps]
nstep = 1000; % number of time steps
A0 = M * 1/dt - A * 1/2; % matrix for step n
A1 = A * 1/2 + M * 1/dt;        % matrix for step n+1
q = qvec/dt;  % source at n=0
phi = A1\q;   % Phi_1
gamma(1,:) = mvec.' * phi; % project to boundary measurements
for i=2:nstep % loop over remaining steps
    q = A0 * phi;
    phi = A1\q;
    gamma(i,:) = mvec.' * phi;
end

plot(dt:dt:dt*nstep,gamma),xlabel('ps'),ylabel('Exitance')