function C = CapacitanceFunction(psi,thetasat,psi0,Az,p)    
%psi: the xylem water potential in Pa
%thetasat= water content at saturation
%psi0= empirical paramter for water pressure of dry xylem
%Az= cross-sectional area for a xylem
%p= empirical curve fitting coefficient for water content
C=Az.*p.*thetasat.*((psi0-psi)./psi0).^(-(p+1));

 