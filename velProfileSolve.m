% This is a workbook just gives velocity profiles

function velProfileSolve
% Fluid Parameters
etaO = 10;
etaI = 0.01;
delta = etaO - etaI;
m = 1;
n = 0.7;

% System Parameters
H = 0.1;
L = 1;
delP = 10;

% Numerical Inputs
itLimit = 10;
DNP = 20;

%Parameter output
parameters = zeros(1, 8);

% Solve for the velocity profile
gammaW = gammaSolve(H, delP, L, etaI, delta, m, n, itLimit);
vP = velPSolve(H, delP, L, etaI, delta, m, n, DNP, itLimit);

parameters(1) = H;
parameters(2) = L;
parameters(3) = m;
parameters(4) = n;
parameters(5) = etaO;
parameters(6) = etaI;
parameters(7) = delP;
parameters(8) = gammaW;

display(vP)
display(parameters)

plotVP(DNP, H, vP)


end


function gamma = gammaSolve(H, delP, L, etaI, delta, m, n, itLimit)
% ////////////////////////////////////////////////////////////////////////
% //////// Here is where we'll do a fixed point iteration
% ////////////////////////////////////////////////////////////////////////

% Numerical Inputs
tol = 1;
iterate = 0;
check = 0;
error = 1e6;

% Initialize
gammaOld = 10;
tau = H*delP/L;

% Carry out fixed point iteration
while (check ~= 1)
    gamma = tau*(etaI + delta/(1 + (m*gammaOld)^n ) )^(-1);
    
    % Checks
    iterate = iterate+1;
    error = abs(gamma - gammaOld);
    if (error < tol) check = 1; end
    if (iterate > itLimit) check = 1; end
    
    gammaOld = gamma;
end
itLimit = 10;
% ////////////////////////////////////////////////////////////////////////
% //////// Here is where finish with Newton's Method
% ////////////////////////////////////////////////////////////////////////

% Numerical Inputs
tol = 1e-10; 
iterate = 0;
check = 0;
error = 1e6;

% Carry out Newton's method
while (check ~= 1)
    f = (etaI + delta/(1 + (m*gammaOld)^n ) )*gammaOld - tau;
    df = etaI + delta/(1 + (m*gammaOld)^n ) ... 
        - n*delta*(m*gammaOld)^n/(1 + (m*gammaOld)^n )^2 ;
    
    gamma = gammaOld - f/df;
    
    % Checks
    iterate = iterate+1;
    error = abs(gamma - gammaOld);
    if (error < tol) check = 1; end
    if (iterate > itLimit) check = 1; end
    
    gammaOld = gamma;
end

end

function v = velSolve(delP, L, etaI, delta, m, n, gammaH, gammaW)
gamma = gammaH;
mGamma = m*gamma;
f = mGamma^n;
g = 1+f;

aF = 1;
bF = 2/n;
cF = (n+1)/n;
dF = -f;
F = real(hypergeom([aF,bF],cF,dF) );

iH = gamma*gamma*(delta*g*F - 2*delta - etaI*g )/(2*g);

gamma = gammaW;
mGamma = m*gamma;
f = mGamma^n;
g = 1+f;

aF = 1;
bF = 2/n;
cF = (n+1)/n;
dF = -f;
F = real(hypergeom([aF,bF],cF,dF) );

iW = gamma*gamma*(delta*g*F - 2*delta - etaI*g )/(2*g);

I = iW - iH;
v = L*I/delP;
end

function vP = velPSolve(H, delP, L, etaI, delta, m, n, DNP, itLimit)
% Lengths
N = DNP;
vLength = 2*DNP -1;

% Calculate values at the wall
gammaW = gammaSolve(H, delP, L, etaI, delta, m, n, itLimit);

velocity = zeros(vLength, 1);
for i = 1:N
    h = (N - i) * H / (N-1);
    gammaH = gammaSolve(h, delP, L, etaI, delta, m, n, itLimit);
    vH = velSolve(delP, L, etaI, delta, m, n, gammaH, gammaW);
    
    velocity(i) = -real(vH);
    if (i < N) velocity(vLength - i + 1) = velocity(i); end
end

vP = velocity;

end

function vPN = velPSolve_Slit_Newtonian(H, L, delP, mu, DNP) 
% Lengths
N = DNP;
vLength = 2*DNP -1;

velocity = zeros(vLength, 1);
for i = 1:N
    h = (N - i) * H / (N-1);
    vH = (H^2 * delP)/(2*mu*L)*(1 - h^2/H^2);
    
    velocity(i) = vH;
    if (i < N) velocity(vLength - i + 1) = velocity(i); end
end

vPN = velocity;
end

function plotVP(DNP, H, vP)
x = zeros(2*DNP-1, 1);
for i=1:DNP
    x(i) = ((i - 1) * H / (DNP-1) );
    x(i+DNP-1) = ((i - 1) * H / (DNP-1) + H ); 
end

v = vP;
plot(x,v);

end