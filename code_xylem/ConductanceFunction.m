function K=ConductanceFunction(psi,Az,Am,ED,kmax,c1,c2)
%psi: the xylem water potential in Pa
%Az:cross-sectional area
%Am:Scaling coefficient
%ED:scaling exponent for effective conductivity (for Norway spruce)
% kmax:maximal xylmem conductance 
%c1:empirical curve fitting coefficient-cavitation pressure
%c2:empirical curve fitting coefficent for conductance
Aaz=(Az/Am).^((ED-2)/2).*Az;   %Effective cross-sectional area
K=Aaz.*kmax.*exp(-(-psi./c1).^c2);

