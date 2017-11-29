function errorCheck
% Fluid inputs
n       =   0.666;
m           =   1;
etaO        =   10;
etaI        =   0.01;
delta       =   etaO-etaI;

% System inputs
H = 1;
L = 10;
delP = 50.05;

% ////////////////////////////////////////////////////////////////////////
% //////// Here is where we'll do a fixed point iteration
% ////////////////////////////////////////////////////////////////////////

% Numerical Inputs
tol = 1;
itLimit = 10; 
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

iterateTotal = iterate;

% ////////////////////////////////////////////////////////////////////////
% //////// Here is where finish with Newton's Method
% ////////////////////////////////////////////////////////////////////////

% Numerical Inputs
tol = 1e-10;
itLimit = 10; 
iterate = 0;
check = 0;
error = 1e6;

% Carry out fixed point iteration
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

iterateTotal = iterateTotal + iterate;

display(error);
display(iterate);
display(iterateTotal);
display(gamma);

end