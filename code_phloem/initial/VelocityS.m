function [ui,vi,pi] = VelocityS(ci,nui,Mu,Psi,X0,G,dz,n)
% pressure
diagP      = - nui(2:n-1).*2/(3*(dz^2)) - Mu;
diagP1     = - Mu; 
diagP2     = - Mu;
% diagP2     = 1;
diagP      = [diagP1; diagP; diagP2];
uppP       = (nui(3:n) - nui(1:n-2))./(12*(dz^2)) + nui(2:n-1)./(3*(dz^2));
uppP       = [0; uppP;0];
lowP       = - (nui(3:n) - nui(1:n-2))./(12*(dz^2)) + nui(2:n-1)./(3*(dz^2));
lowP       = [0;lowP; 0];
eqP        = - ci(2:n-1) - X0.*Psi(2:n-1) + G*dz.*((2:n-1)' - 1/2);
eqP1       = - ci(1)  -  X0*Psi(1) + G*dz*0.5;
eqP2       = - ci(n) - X0*Psi(n) + G*dz*(n - 1/2);
% eqP2       = 0;
eqP        = [eqP1;eqP;eqP2];
pi         = Thomas(lowP,diagP,uppP,eqP);
pi         = pi';
% axial velocity
ui         = (nui(2:n) + nui(1:n-1)).*( pi(1:n-1) - pi(2:n) )./(6*dz);
% radial velocity
vi         = zeros(n,1);
vi(2:n-1)  = (3/8)*dz.*(ui(2:n-1) - ui(1:n-2));
%%%%%%%%%% BC P = 0
% vi(n)      = - (3/8)*( G*dz*(n-1/2) - ci(n) - X0*Psi(n) );
%%%%%%%%%%
end