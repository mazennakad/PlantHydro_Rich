function [NP,NU,NV,NC,NN,F5] = BuildN(n,cw,bb,c,nu)
% P variable
NP = zeros(n);
% U variable
NU = zeros(n,n-1);
% V variable
NV = zeros(n);
% C variable 
n1 =  - ( - bb(2)*cw - 2*bb(3)*(cw^2).*c ...
    - 3*bb(4)*(cw^3).*c.^2 - 4*bb(5)*(cw^4).*c.^3 ) ...
   .*exp( - bb(2)*cw.*(c - 1) - bb(3)*(cw^2).*(c.^2 - 1) ...
    - bb(4)*(cw^3).*(c.^3 - 1) - bb(5)*(cw^4).*(c.^4 - 1) );
NC = diag(n1);
% Nu variable
NN = diag(ones(1,n));

% Equality
F5 = nu - exp( - bb(2)*cw.*(c - 1) - bb(3)*(cw^2).*(c.^2 - 1) ...
    - bb(4)*(cw^3).*(c.^3 - 1) - bb(5)*(cw^4).*(c.^4 - 1) );
end