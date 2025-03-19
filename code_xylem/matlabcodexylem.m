n=5;  %Number of grids
L=2;   %length of the xylem
g=9.81;  %constant of gravitational acceleration
rho=1000;  %water density
dz=L/n;  % grid size
dt=0.05;
t=0:dt:200*dt; % duration of simulation
T = 0.1; % transpiration

% capacitance constants
thetasat=573.5;  %water content at saturation
psi0=5.74*10^8;  %empirical paramter for water pressure of dry xylem
R=0.5;      %typical radius for a xylem
Az=pi*R.^2; %cross-sectional area for a xylem
p=20;       %empirical curve fitting coefficient for water content

% Conductance constants
Am = 0.01;
ED = 2.44;
kmax=1.36*10^(-8); %maximal xylmem conductance 
c1=4.8^10^(-6); %empirical curve fitting coefficient-cavitation pressure
c2=3.5;   %empirical curve fitting coefficent for conductance



psis=-783.77;   %water pressure at the base of the xylem




psik = zeros(1,n);
psik(1) = psis;
psi = psik;


C = Fcapacitance(psi,thetasat,psi0,Az,p);
K = Fconductance(psi,Az,Am,ED,kmax,c1,c2);
F_values=F(psi,psik,n,g,rho,dz,dt,C,K,T);
%dpsi = JFNK();
%for i = 1:200
%    while check == true
%    end
%end



