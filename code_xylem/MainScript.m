n=1000;  %Number of grids
L=2;   %length of the xylem(m)
g=9.81;  %constant of gravitational acceleration(m/s^2)
rho=1000;  %water density(kg/m^3)
dz=L/n;  % grid size
dt=0.05; %time step(s)
t=0:dt:200*dt; % duration of simulation
T =0.01; % transpiration(kg/s)

% capacitance constants
thetasat=573.5;  %water content at saturation(kg/m^3)
psi0=5.74*10^8;  %empirical paramter for water pressure of dry xylem(Pa)
R=0.01;      %typical radius for a xylem(m)
Az=pi*R.^2; %cross-sectional area for a xylem
p=20;       %empirical curve fitting coefficient for water content

% Conductance constants
Am = 0.01;   %scaling coefficient(m^2)
ED = 2.44;   %exponent order coefficient
kmax=1.36*10^(-8); %maximal xylmem conductance (m^2/s)
c1=4.8*10^(-6); %empirical curve fitting coefficient-cavitation pressure(Pa)
c2=3.5;   %empirical curve fitting coefficent for conductance
n=5;
psis=-783.77;   %water pressure at the base of the xylem (Pa)
%psik = psis.*ones(n,1);   %psi at iteration 1 at the base of the xylem
psik = linspace(psis,0,n);
%psik(1) = psis;
psi = psik;
C = CapacitanceFunction(psi,thetasat,psi0,Az,p);
K = ConductanceFunction(psi,Az,Am,ED,kmax,c1,c2);
F_values=MainFunction(psi,psik,n,g,rho,dz,dt,C,K,T);
%dpsi = JFNK();
%initializing parameters


m=500;
h=zeros(m+1,m);  %Hessenberg matrix
w=zeros(n,1);  %w=Jv
Vm=zeros(n,m);
eps = 1e-8;
e1=zeros(m+1,1);
e1(1)=1;
r = -(MainFunction(psi,psik,n,g,rho,dz,dt,C,K,T));
norm(r)
while norm(r)>1e-9
for i=1:m
    if i==1
        %r=-(MainFunction(psi,psik,n,g,rho,dz,dt,C,K,T)); %residual
        rn=norm(r);
        v = r./rn;  
    end
    Vm(:,i) = v;
    w=(MainFunction(psi+eps.*v,psik,n,g,rho,dz,dt,C,K,T)-MainFunction(psi,psik,n,g,rho,dz,dt,C,K,T))./eps;
    wp = w;

    %h(i,i)=Vm(:,i)'.*w;
    %wp= w - h(i,i).*Vm(:,i);
    %h(i+1,i) = norm(wp);
    %v = wp./h(i+1,i);

    for j=1:i
        h(j, i) = Vm(:, j)' * w;
        wp=wp-h(j,i).*Vm(:,j);
    end
    h(i+1,i)=norm(wp);
    v = wp./h(i+1,i);
    %if h(i + 1, i) ~= 0
    %    v = wp ./ h(i + 1, i);
    %else
    %    break
    %end
end
y = h \ (-rn.*e1(1:i+1));
% y=-rn.*e1(1:i+1)./h(1:i+1,1:i);
deltap=mtimes(Vm,y);
psi=psi+deltap;
C = CapacitanceFunction(psi,thetasat,psi0,Az,p);
K = ConductanceFunction(psi,Az,Am,ED,kmax,c1,c2);
r=-(MainFunction(psi,psik,n,g,rho,dz,dt,C,K,T));
norm(r)
end




% check=true;
% 
% 
%     v(1)=r/norm(r);
%     for m=1:5
%     w=(F(psi+dt.*v,psik,n,g,rho,dz,dr,C,K,T)-F(psi,psik,n,g,rho,dz,dt,C,K,T))/dt;
%     h(m,m)=v(m).*w;
%     for j=1:m
%     wp=w;
%     wp=wp-h(j,m).*v(j);
%     h(m+1,m)=norm(wp);
%     end;
%     v(m)=wp/h(m+1,m);
%     end;
% 
% 
% 
% for i = 1:200
%     while check== true
%     r=-F(psi,pisk,thetasat,psi0,Az,p);
%     v=zeros(1,5)
%     v(1)=r/norm(r);
%     for m=1:5
%     w=(F(psi+dt.*v,psik,n,g,rho,dz,dr,C,K,T)-F(psi,psik,n,g,rho,dz,dt,C,K,T))/dt;
%     h(m,m)=v(m).*w;
%     for j=1:m
%     wp=w;
%     wp=wp-h(j,m).*v(j);
%     h(m+1,m)=norm(wp);
%     end;
%     v(m)=wp/h(m+1,m);
%     end;
% 
% 
% v(i+1)=wp/(h(i+1,i));
% if norm(r)<10^(-10)
%     check==false;
% %    end
% %end



