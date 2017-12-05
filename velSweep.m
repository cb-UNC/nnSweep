%{
    This code carries out a sweep of velocity profiles across a set of
    fluids with differing Cross model parameters. The velocity profile
    flattens out as delta gets larger, where delta = etaO - etaI

    The form of the cross model used here is:

    eta = etaI - delta/(1 + (m*gamma)^n);

    where
    eta = viscosity
    etaI = infinite shear viscosity
    etaO = zero shear viscosity
    delta = etaO - etaI
    m = molecular constant
    n = behavior index
    gamma = shear rate

    There is a starting eta, a maximum etaO, and a minimum etaI, calculated
    by the following equations:

    N = sweepNumber-1;
    etaOScale = (log10(etaOMax) - log10(etaStart))/N;
    etaIScale = (log10(etaStart) - log10(etaIMin))/N;
    etaO = etaStart * 10^(etaOScale * i);
    etaI = etaStart * 10^(- etaIScale * i);

    These eta parameters are calculated for each sweep, m, n are constant.

    For the slit, the H, L, and gamma at the wall are held constant, the
    change in pressure is calculated by calculating the eta at the wall,
    then:

    delP = etaW * gammaW * L / H;

    This code outputs a vector, sweepVPs, which holds all of the velocity
    profiles, and it plots the velocity profiles normalized by their
    largest velocity. 

    See code below for parameter values.
%}

function velSweep
% Eta scaling parameters
etaStart = 1.0;
etaOMax = 10;
etaIMin = 0.01;
sweepNumber = 8;

N = sweepNumber-1;
etaOScale = (log10(etaOMax) - log10(etaStart))/N;
etaIScale = (log10(etaStart) - log10(etaIMin))/N;

% Constant parameters
m = 100.0;
n = 0.7;
H = 1;
L = 10;
gammaMax = 10^2;
itLimit = 10;
DNP = 9;

sweepVPs = zeros(2*DNP-1, sweepNumber);
for i = 1:sweepNumber
    % Adjusted parameters
    etaO = etaStart * 10^(etaOScale * i); 
    etaI = etaStart * 10^(- etaIScale * i);
    delta = etaO - etaI; 
    gammaW = gammaMax; 
    etaW = etaI + delta/(1 + (m*gammaW)^n); 
    delP = etaW * gammaW * L / H;

    vP = velPSolve(H, delP, L, etaI, delta, m, n, DNP, itLimit);
    for j = 1:DNP
        sweepVPs(j,i) = vP(j);
        if (i < DNP) sweepVPs(2*DNP - j, i) = vP(j); end
    end
end

sweepVPs
x = zeros(2*DNP-1, 1);
for i=1:DNP
    x(i) = (i - 1) * H / (DNP-1);
    x(i+DNP-1) = (i - 1) * H / (DNP-1) + H; 
end

hold on
v = zeros(2*DNP-1,1);
for i = 1:sweepNumber
    vMax = sweepVPs(DNP,i);
    for j = 1:(2*DNP-1)
        v(j) = sweepVPs(j,i)/vMax;
    end
    plot(x,v)
end

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