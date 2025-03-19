function [F_values]=F(psi,psip,n,g,rho,dz,dt,C,K,T)

% psi is the xylem water potential is Pa
% psip is the xylem water potential at previous iteration
% L is the stem length
% n is the number of grids
% g is the gravitational constant
% dz is the grid size
% dt is the time step

F_values = zeros(1,n);
F_values(n) = C(n).*( psi(n) - psip(n) )./dt ...
    - T./dz + K(n-1).*( ( psi(n)-psi(n-1) )./dz +rho*g )./dz;
F_values(2:n-1) = C(2:n-1).*( psi(2:n-1) - psip(2:n-1) )./dt ...
    - K(2:n-1).*( ( psi(3:n)-psi(2:n-1) )./dz + rho.*g )./dz ...
    + K(1:n-2).*( ( psi(2:n-1)-psi(1:n-2) )./dz +rho*g )./dz;

