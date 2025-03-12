function [ci,ui,vi,pi,nui] = guess(n,dz,G,Mu,X0,Psi,cw,bb,m)
% concentration
co         = 2.*(G*dz.*( (1:n)' - 1/2 ) - X0.*Psi);

switch m
    case 0
        nuo = ones(n,1);
    case 1
        nuo = exp( - bb(2)*cw.*(co - 1) ...
            - bb(3)*(cw^2).*(co.^2 - 1) ...
            - bb(4)*(cw^3).*(co.^3 - 1) ...
            - bb(5)*(cw^4).*(co.^4 - 1) );
end
[uo,vo,po] = VelocityS(co,nuo,Mu,Psi,X0,G,dz,n);
ci = co;
ui = uo;
vi = vo;
pi = po;
nui = nuo;

end