function C1 = PriorOU(mesh,std_,r1)

g1 = mesh.g;
N1 = size(g1,1);
C1 = zeros(N1,N1); % the covariance matrix

for i=1:N1
    Rsq = sqrt((g1(:,1) - g1(i,1)).^2 + (g1(:,2) - g1(i,2)).^2); 
    C1(i,:) = std_^2*exp(-Rsq/r1);   
end
