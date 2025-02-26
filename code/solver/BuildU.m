function [UP,UU,UV,UC,UN,F2] = BuildU(n,dz,u,p,nu)

% P variable
UP = diag( - (nu(2:n) + nu(1:n-1))./(6*dz) ) ...
    + diag( (nu(2:n-1) + nu(1:n-2))./(6*dz),1 );
UP = [UP zeros(n-1,1)];
UP(n-1,n) = (nu(n) + nu(n-1))/(6*dz);

% U variable
UU = diag(ones(1,n-1)); 
% V variable 
UV = zeros(n-1,n);
% C variable
UC = zeros(n-1,n);
% Nu variable
UN = diag( (p(2:n) - p(1:n-1))./(6*dz) ) ...
    + diag( (p(2:n-1) - p(1:n-2))./(6*dz),1 );
UN = [UN zeros(n-1,1)];
UN(n-1,n) = (p(n) - p(n-1))/(6*dz);

% Equality
F2 = u + (nu(2:n) + nu(1:n-1)).*(p(2:n) - p(1:n-1))./(6*dz);
end