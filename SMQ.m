clear all;
clc;
n = 10000; % number of nodes
nsolve = 20;
x0 = -50;
x1 = 50;
h = (x1-x0)/(n-1);
x = linspace(x0,x1,n);

e = (1/h^2)*ones(n-2,1); % n-2 is used to pull off the boundary conditions
A = spdiags([-e 2*e -e],[-1 0 1],n-2,n-2); % central difference corresponding to -0.5*laplacian

x = x(2:end-1); % pull off the boundary conditions
c = 0.5;
V = c*(x.^2-x).*(abs(x-0.5)<=0.5); % harmonic potential
bottom = min(V);
%[x' V']
%V = 0.5*x.^2;
V = spdiags(V',0,n-2,n-2);
H = A+V;
[U E] = eigs(H,nsolve,'sm');
plot(x,diag(V),'k-');
hold on;
for i = 0:nsolve-1
    plot(x,-U(:,nsolve-i)+E(nsolve-i,nsolve-i));
end
diag(E)
