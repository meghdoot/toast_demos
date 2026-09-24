% Demo 1: data generation, homogeneous propert

clc, clear all, close all,

% create mesh
[vtx,idx,eltp] = mkcircle(25,6,32,2); % radius, sectors, rings, bdd
% create mesh object
mesh = toastMesh (vtx,idx,eltp);
% specify sources and detectors
ths = (-pi:2*pi/8:pi)';thd = 0.4-ths; 
sp = [25*cos(ths) 25*sin(ths)];mp = [25*cos(thd) 25*sin(thd)];

figure,
mesh.Display(),hold on,
plot(sp(:,1),sp(:,2),'r*'),plot(mp(:,1),mp(:,2),'b*')

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
Phi = K\qvec; % Fluence
Gamma = mvec.' * Phi; % Exitance

% Display Phi and Gamma
figure,
subplot(121),mesh.Display(log(abs(Phi(:,1)))),title('LAmp')% log-amplitude
subplot(122),mesh.Display(angle(Phi(:,1))),title('Phs')% phase

figure,
plot([log(abs(Gamma(:)));angle(Gamma(:))]),title('Data')




