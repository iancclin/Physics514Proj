clear all;
clc;
m = 1;
hbar = 1;
c = 0.5;
V = @(x) c*(x.^2-x).*(abs(x-0.5)<=0.5);%(0.5*x.^2).*(abs(x)<=5);

%%  define domain
xR = 1; xL = 0;
nodes = 10000;
x = linspace(xL,xR,nodes);
h = (xR-xL)/(nodes-1);
%b = 0.23;
%nodeb = round((b-xL)/(xR-xL)*nodes)+1;
%b = x(nodeb);

test = zeros(length(x),1);
psi = zeros(length(x)-2,length(x));

for i = 2:length(x)-1
    b = x(i);
    
    pot = V(x);
    E = V(b);
    k = (2*m/(hbar^2))*(pot-E);
    f = (h^2/12)*k;
    fL0 = (h^2/12)*(2*m/(hbar^2))*(V(xL-h)-E);
    fR0 = (h^2/12)*(2*m/(hbar^2))*(V(xR+h)-E);
    
    
    %% define boundary function values
    psiR(end) = exp(-0.5);
    psiR0 = exp(-0.5-h);
    psiL(1) = exp(-0.5);
    psiL0 = exp(-0.5-h);
    
    
    
    %% calculate the first 2 values
    psiL(2) = (1/(1-f(2)))*(2*(1+5*f(1))*psiL(1)-(1-fL0)*psiL0);
    psiR(end-1) = (1/(1-f(end-1)))*(2*(1+5*f(end))*psiR(end)-(1-fR0)*psiR0);
    for n = 3:i
        psiL(n) = (1/(1-f(n)))*(2*(1+5*f(n-1))*psiL(n-1)-(1-f(n-2))*psiL(n-2));
    end
    for n = nodes-2:-1:i
        psiR(n) = (1/(1-f(n)))*(2*(1+5*f(n+1))*psiR(n+1)-(1-f(n+2))*psiR(n+2));
    end
    diffPsibL = (psiL(i)-psiL(i-1))/h;
    diffPsibR = (psiR(i+1)-psiR(i))/h;
    test(i) = diffPsibR/psiR(i) - diffPsibL/psiL(i);
    psi(i,:) = [psiL(1:i)' psiR(i+1:end)'];
end
res = [test x'];
