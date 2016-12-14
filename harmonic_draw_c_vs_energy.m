clear all;
clc;

%% define quantum constants
m = 1;
hbar = 1;

%%  define domain and meshe size
% for harmonic test use
% xL = -5; xR = 5;
xL = 0; xR = 1;
nodes = 10000;

% define how many states we want to investigate
states = 0:1:5;
% define the constant for the potential well
c = 0:1:400;

for s = 1:length(states)
    state = states(s);
    eplot = zeros(1,length(c));
    
    for k = 1:length(c)
        V = @(x) c(k)*(x.^2-x).*(abs(x-0.5)<=0.5);
        % harmonic well used for testing Numerov algorithm
        % V = @(x) (0.5*x.^2);
        
        %% define the start points and allocate the variable
        % start from finding the ground state
        psi = [];
        energy = [];
        position = [];
        % The gournd staet will be, in general, at neither the top nor the bottom
        % of the potential well. As a result, initialize it as the point in
        % between.
        
        xmin = 0.5*(xL+xR);
        Emin = V(xmin);
        xmax = xL;
        Emax = V(xmax);
        b = 0.5*(xmin+xmax);
        
        tol = 1e-6;
        % for bisection search, it is impossible to search more than log2(nodes)
        % times, otherwise memory violation will happen
        maxIter = 100;
        
        iter = 1;
        toldiff = 1;
        
        while iter <= maxIter
            iter = iter + 1;
            E = V(b);
            nodesL = round((b-xL)/(xR-xL)*nodes)+1;
            nodesR = round((xR-b)/(xR-xL)*nodes);
            if nodesL <=2 || nodesR <= 2
                fprintf('unable to find the solution of state %d after %d iterations.\n', state, iter);
                break;
            end
            gridL = linspace(xL,b,nodesL);
            hL = gridL(2)-gridL(1);
            potL = V(gridL);
            kL = (2*m/(hbar^2))*(potL-E);
            fL = (hL^2/12)*kL;
            fL0 = (hL^2/12)*(2*m/(hbar^2))*(V(xL-hL)-E);
            psiL = zeros(1,nodesL);
            
            gridR = linspace(b,xR,nodesR);
            hR = gridR(2)-gridR(1);
            potR = V(gridR);
            kR = (2*m/(hbar^2))*(potR-E);
            fR = (hR^2/12)*kR;
            fR0 = (hR^2/12)*(2*m/(hbar^2))*(V(xR+hR)-E);
            psiR = zeros(1,nodesR);
            
            %% define boundary function values
            psiR(end) = exp(-0.5);
            psiR0 = exp(-0.5-hR);
            psiL(1) = exp(-0.5);
            psiL0 = exp(-0.5-hL);
            
            %% calculate the first 2 values
            psiL(2) = (1/(1-fL(2)))*(2*(1+5*fL(1))*psiL(1)-(1-fL0)*psiL0);
            psiR(end-1) = (1/(1-fR(end-1)))*(2*(1+5*fR(end))*psiR(end)-(1-fR0)*psiR0);
            
            signChange = 0;
            
            for n = 3:nodesL
                psiL(n) = (1/(1-fL(n)))*(2*(1+5*fL(n-1))*psiL(n-1)-(1-fL(n-2))*psiL(n-2));
                if psiL(n)*psiL(n-1) < 0
                    signChange = signChange + 1;
                end
            end
            for n = nodesR-2:-1:1
                psiR(n) = (1/(1-fR(n)))*(2*(1+5*fR(n+1))*psiR(n+1)-(1-fR(n+2))*psiR(n+2));
                if psiR(n)*psiR(n+1) < 0
                    signChange = signChange + 1;
                end
            end
            % too many nodes, meaning the energy is too high
            if signChange > state
                Emax = E;
                xmax = b;
                b = 0.5*(xmin+xmax);
            elseif signChange < state
                Emin = E;
                xmin = b;
                b = 0.5*(xmin+xmax);
            else
                diffPsibL = (psiL(end)-psiL(end-1))/hL;
                diffPsibR = (psiR(2)-psiR(1))/hR;
                % if diffVal < 0, then the guess energy is too high and vice versa
                diffVal = diffPsibR/psiR(1) - diffPsibL/psiL(end);
                if abs(diffVal) < tol
                    energy(state+1) = E;
                    position(state+1,:) = [gridL gridR(2:end)];
                    psiR = psiR*(psiL(end)/psiR(1));
                    intL = hL*sum((0.5*(psiL(1:end-1)+psiL(2:end))).^2);
                    intR = hR*sum((0.5*(psiR(1:end-1)+psiR(2:end))).^2);
                    normPsi = intL+intR;
                    psi(state+1,:) = [psiL psiR(2:end)]/normPsi;
                    break;
                elseif diffVal > 0
                    Emax = E;
                    xmax = b;
                    b = 0.5*(xmin+xmax);
                else
                    Emin = E;
                    xmin = b;
                    b = 0.5*(xmin+xmax);
                end
            end
            
            %res = [test x'];
        end
        %plot(position(state+1,:),psi(state+1,:));
        if iter > maxIter
            fprintf('unable to find the solution of state %d.\n', state);
        end
        if isempty(energy)
            eplot(k) = 0;
        else
            eplot(k) = energy(state+1);
        end
    end
    plot(c,eplot);
    hold on;
end
legend('ground state','1st excited state','2nd excited state','3rd excited state');
