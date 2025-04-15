function [F_values]=MainFunction(psi,psip,n,g,rho,dz,dt,C,K,T)

% psi: the xylem water potential in Pa
% psip: xylem water potential at previous time
% L: stem length
% n: number of grids
% g:the gravitational constant
% dz:the grid size
% dt:the time step

F_values = zeros(n,1);
F_values(n) = C(n).*( psi(n) - psip(n) )./dt ...
    - T./dz + K(n-1).*( ( psi(n)-psi(n-1) )./dz +rho*g )./dz;
F_values(2:n-1) = C(2:n-1).*( psi(2:n-1) - psip(2:n-1) )./dt ...
    - K(2:n-1).*( ( psi(3:n)-psi(2:n-1) )./dz + rho.*g )./dz ...
    + K(1:n-2).*( ( psi(2:n-1)-psi(1:n-2) )./dz +rho*g )./dz;