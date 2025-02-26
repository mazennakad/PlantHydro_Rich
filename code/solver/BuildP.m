function [PP,PU,PV,PC,PN,F1] = BuildP(n,dz,Mu,G,X0,Psi,c,p,nu)

% P variable
diagP  = - nu(2:n-1).*2/(3*(dz^2)) - Mu;
diagP1 = - Mu; 
diagP2 = - Mu;
% diagP2 = 1;
diagP  = [diagP1; diagP; diagP2];
uppP   = (nu(3:n) - nu(1:n-2))./(12*(dz^2)) + nu(2:n-1)./(3*(dz^2));
uppP   = [0; uppP];
lowP   = - (nu(3:n) - nu(1:n-2))./(12*(dz^2)) + nu(2:n-1)./(3*(dz^2));
lowP   = [lowP; 0];
PP     = diag(diagP) + diag(uppP,1) + diag(lowP,-1);
% U variable
PU     = zeros(n,n-1);
% V variable
PV     = zeros(n);
% C variable
PC     = diag(ones(1,n));
%%%%%%%% P = 0
% PC(end)= 0; 
%%%%%%%%
% Nu variable
diagN  = ( p(3:n) - 2.*p(2:n-1) + p(1:n-2) )./(3*(dz^2));
diagN  = [0; diagN; 0];
uppN   = (p(3:n) - p(1:n-2))./(12*(dz^2));
lowN   = - uppN;
uppN   = [0; uppN];
lowN   = [lowN; 0];
PN     = diag(diagN) + diag(uppN,1) + diag(lowN,-1);
% Equality
F1     = (nu(3:n) - nu(1:n-2)).*(p(3:n) - p(1:n-2))./(6*(dz^2)) ...
       + nu(2:n-1).*( p(3:n) - 2.*p(2:n-1) + p(1:n-2) )./(3*(dz^2)) ...
       - Mu.*p(2:n-1) + c(2:n-1) + X0.*Psi(2:n-1) - G*dz.*((2:n-1)' - 1/2);
F11    = - Mu*p(1) + c(1) + X0*Psi(1) - G*dz*0.5;
F12    = - Mu*p(n) + c(n) + X0*Psi(n) - G*dz*(n - 1/2);
% F12    = 0;
F1     = [F11; F1; F12];
end