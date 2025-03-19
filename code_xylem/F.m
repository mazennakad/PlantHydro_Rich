function [F_values]=F(psi,psizero)
L=2;   %length of the xylem
i=5;   %number of sections into which the tree is divided
g=9.81;  %constant of gravitational acceleration
po=999;  %water density
dz=L/i;  %distance between each iteration
dt=0.05;
F_values=zeros(1,i)
for j=2:1:i
if j==1
    psi(j-1)=0;
elseif j==i
    psi(j+1)=0
end
C=Fcapacitance(psi(j))
K=Fconductance(psi(j))
F_values(j)=C.*((psi(j)-psizero(j)))/dt-((K.*(psi(j+1)-psi(j))./dz+po.*g)-(K.*(psi(j)-psi(j-1))./dz+po.*g))./dz
end