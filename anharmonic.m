clc;
n = 6;
arr = sqrt(1:n)';
aan = spdiags(sqrt(0:n-1)',1,n,n);
acr = spdiags(sqrt(1:n)',-1,n,n);
x = (1/sqrt(2))*(acr+aan);
syms lam;
H0 = (acr*aan)+spdiags(0.5*ones(n,1),0,n,n);

lam=0:0.01:1;
numeig = 20;
val = zeros(numeig,length(lam));
valex = val;

for i = 1:length(lam)
    H = H0+lam(i)*x^4;
    val(:,i) = eigs(H,numeig,'SM');
    ex = eig(H);
    valex(:,i) = ex(1:numeig);
end
valex(1:end,:) = valex(end:-1:1,:);
plot(lam,valex(1,:),'k-');
hold on;
plot(lam,val(1,:),'b--');
