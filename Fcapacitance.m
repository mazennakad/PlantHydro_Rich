function C=Fcapacitance(psi)
thetasat=573.5;  %water content at saturation
psi0=5.74*10^8;  %empirical paramter for water pressure of dry xylem
R=0.5;  %typical radius for a xylem
Az=pi*R.^2; %cross-sectional area for a xylem
p=20;  %empirical curve fitting coefficient for water content
C=Az.*p.*thetasat.*((psi0-psi)./psi0).^(-(p+1));

 