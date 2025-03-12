function [VP,VU,VV,VC,VN,F3] = BuildV(n,dz,v,p,nu,c,G,X0,Psi)
% P variable
diagP = - nu(2:n-1)./(4*(dz^2));
uppP  = (nu(3:n) - nu(1:n-2))./(32*(dz^2)) ...
      + nu(2:n-1)./(8*(dz^2));
lowP  = - (nu(3:n) - nu(1:n-2))./(32*(dz^2)) ...
      + nu(2:n-1)./(8*(dz^2));
diagP = [0; diagP; 0];
uppP  = [0; uppP];
lowP  = [lowP; 0];
VP    = diag(diagP) + diag(lowP,-1) + diag(uppP,1);
% U variable
VU    = zeros(n,n-1); 
% V variable 
VV    = diag(ones(1,n)); 
% C variable
VC    = zeros(n);
% Nu variable
diagN = (p(3:n) - 2.*p(2:n-1) + p(1:n-2))./(8*(dz^2));
uppN  = (p(3:n) - p(1:n-2))./(32*(dz^2));
lowN  = - uppN;
diagN = [0; diagN; 0];
uppN  = [0; uppN];
lowN  = [lowN; 0];
VN    = diag(diagN) + diag(lowN,-1) + diag(uppN,1);
% Equality
F3 = v(2:n-1) ...
    + (nu(3:n) - nu(1:n-2)).*(p(3:n) - p(1:n-2))./(32*(dz^2))...
    + nu(2:n-1).*(p(3:n) - 2.*p(2:n-1) + p(1:n-2))./(8*(dz^2));
F3 = [0; F3; 0];
% F3 = [0; F3; v(n) + (3/8)*(G*dz*(n-1/2) - c(n) - X0*Psi(n))];
end