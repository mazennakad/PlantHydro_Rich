function [CP,CU,CV,CC,CN,F4] = BuildC(n,dz,Pe,Sl,Ss,Sr,es,c,u,v,s,dt,cp,k1,Ss1,Ss2,nts)
% P variable       
CP = zeros(n); 

switch k1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Constant sink %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 0
        % U variable
        a1        = ( Pe - (7/45)*(Pe^2).*v(1:n-1) ).*c(1:n-1)./dz ...
                  - (4/105)*((Pe/dz)^2).*u.*(c(2:n) - c(1:n-1)) ...
                  - (7/120)*(Ss/es)*(Pe^2)/dz;
        a2        = - ( Pe - (7/45)*(Pe^2).*v(1:n-2) ).*c(1:n-2)./dz ...
                  + (4/105)*((Pe/dz)^2).*u(1:n-2).*(c(2:n-1) - c(1:n-2)) ...
                  + (7/120)*(Ss/es)*(Pe^2)/dz;
        CU        = diag(a1) + diag(a2,-1);
        CU        = [CU; zeros(1,n-1)];
        CU(n,n-1) =  - ( Pe - (7/45)*(Pe^2).*v(n-1) ).*c(n-1)/dz ...
                  + (4/105)*((Pe/dz)^2).*u(n-1).*(c(n) - c(n-1)) ...
                  + (7/120)*(Ss/es)*(Pe^2)/dz;
        % V variable
        a3        = - (7/45)*(Pe^2).*u(1:n-1).*c(1:n-1)./dz;
        a4        = (7/45)*(Pe^2).*u(1:n-2).*c(1:n-2)./dz;
        CV        = diag(a3) + diag(a4,-1);
        CV        = [CV; zeros(1,n-1)];
        CV        = [CV, zeros(n,1)];
        CV(n,n-1) = (7/45)*(Pe^2).*u(n-1).*c(n-1)./dz;
        CV(1,1)   = 0;

        switch s
            case 1
                 % C variable
                 diagC  = 1/dt + ( Pe - (7/45)*(Pe^2).*v(2:n-1) ).*u(2:n-1)./dz ...
                        + ( (2/105)*(Pe^2).*( u(2:n-1).^2 ) + es^2 )./(dz^2)...
                        + ( (2/105)*(Pe^2).*( u(1:n-2).^2 ) + es^2 )./(dz^2);
                 diagC1 = 1/dt + Pe*u(1)/dz ...
                        + ( (2/105)*(Pe^2).*( u(1).^2 ) + es^2 )./(dz^2);
                 diagC2 = 1/dt + ( (2/105)*(Pe^2).*( u(n-1).^2 ) + es^2 )./(dz^2);
                 diagC  = [diagC1; diagC; diagC2];
                 uppC   = - ( (2/105)*(Pe^2).*( u(1:n-1).^2 ) + es^2 )./(dz^2);
                 lowC   = - ( Pe - (7/45)*(Pe^2).*v(1:n-1) ).*u./dz ...
                        - ( (2/105)*(Pe^2).*( u.^2 ) + es^2 )./(dz^2);
                 CC     = diag(diagC) + diag(uppC,1) + diag(lowC,-1);
                 % Nu variable
                 CN     = zeros(n);
                 % Equality
                 F4     = (c(2:n-1) - cp(2:n-1))./dt ...
                        + (Pe - (7/45)*(Pe^2).*v(2:n-1)).*u(2:n-1).*c(2:n-1)./dz...
                        - ( (2/105)*(Pe^2).*(u(2:n-1).^2) + es^2 ).*(c(3:n) - c(2:n-1))./(dz^2)...
                        + ( (Pe*dz/es).*(2:n-1)' - (7/120)*(Pe^2).*u(2:n-1)./es ).*Ss./dz...
                        - (Pe - (7/45)*(Pe^2).*v(1:n-2)).*u(1:n-2).*c(1:n-2)./dz...
                        + ( (2/105)*(Pe^2).*(u(1:n-2).^2) + es^2 ).*(c(2:n-1) - c(1:n-2))./(dz^2)...
                        - ( (Pe*dz/es).*(1:n-2)' - (7/120)*(Pe^2).*u(1:n-2)./es ).*Ss./dz;
                 F41    = (c(1) - cp(1))./dt ...
                        + (Pe - (7/45)*(Pe^2).*v(1)).*u(1).*c(1)./dz...
                        - ( (2/105)*(Pe^2).*(u(1).^2) + es^2 ).*(c(2) - c(1))./(dz^2)...
                        + ( (Pe*dz/es) - (7/120)*(Pe^2).*u(1)./es ).*Ss./dz...
                        - Pe*Sl/dz;
                 F42    = (c(n) - cp(n))./dt ...
                        + (Pe*Sr + (Pe/es)*Ss)/dz ...
                        - (Pe - (7/45)*(Pe^2).*v(n-1)).*u(n-1).*c(n-1)./dz...
                        + ( (2/105)*(Pe^2).*(u(n-1).^2) + es^2 ).*(c(n) - c(n-1))./(dz^2)...
                        - ( (Pe*dz/es).*(n-1) - (7/120)*(Pe^2).*u(n-1)./es ).*Ss./dz;
                 F4     = [F41; F4; F42];
                 

            case 0
                % C variable
                 diagC  = ( Pe - (7/45)*(Pe^2).*v(2:n-1) ).*u(2:n-1)./dz ...
                        + ( (2/105)*(Pe^2).*( u(2:n-1).^2 ) + es^2 )./(dz^2)...
                        + ( (2/105)*(Pe^2).*( u(1:n-2).^2 ) + es^2 )./(dz^2);
                 diagC1 = Pe*u(1)/dz ...
                        + ( (2/105)*(Pe^2).*( u(1).^2 ) + es^2 )./(dz^2);
                 diagC2 = ( (2/105)*(Pe^2).*( u(n-1).^2 ) + es^2 )./(dz^2);
                 % diagC2 = 1;
                 diagC  = [diagC1; diagC; diagC2];
                 uppC   = - ( (2/105)*(Pe^2).*( u(1:n-1).^2 ) + es^2 )./(dz^2);
                 lowC   = - ( Pe - (7/45)*(Pe^2).*v(1:n-1) ).*u./dz ...
                        - ( (2/105)*(Pe^2).*( u.^2 ) + es^2 )./(dz^2);
                 CC     = diag(diagC) + diag(uppC,1) + diag(lowC,-1);
                 % Nu variable
                 CN     = zeros(n);
                 % Equality
                 F4     = (Pe - (7/45)*(Pe^2).*v(2:n-1)).*u(2:n-1).*c(2:n-1)./dz...
                        - ( (2/105)*(Pe^2).*(u(2:n-1).^2) + es^2 ).*(c(3:n) - c(2:n-1))./(dz^2)...
                        + ( (Pe*dz/es).*(2:n-1)' - (7/120)*(Pe^2).*u(2:n-1)./es ).*Ss./dz...
                        - (Pe - (7/45)*(Pe^2).*v(1:n-2)).*u(1:n-2).*c(1:n-2)./dz...
                        + ( (2/105)*(Pe^2).*(u(1:n-2).^2) + es^2 ).*(c(2:n-1) - c(1:n-2))./(dz^2)...
                        - ( (Pe*dz/es).*(1:n-2)' - (7/120)*(Pe^2).*u(1:n-2)./es ).*Ss./dz;
                 F41    = Pe.*u(1).*c(1)./dz...
                        - ( (2/105)*(Pe^2).*(u(1).^2) + es^2 ).*(c(2) - c(1))./(dz^2)...
                        + ( (Pe*dz/es) - (7/120)*(Pe^2).*u(1)./es ).*Ss./dz...
                        - Pe*Sl/dz;
                 F42    = (Pe*Sr + (Pe/es)*Ss)/dz ...
                        - (Pe - (7/45)*(Pe^2).*v(n-1)).*u(n-1).*c(n-1)./dz...
                        + ( (2/105)*(Pe^2).*(u(n-1).^2) + es^2 ).*(c(n) - c(n-1))./(dz^2)...
                        - ( (Pe*dz/es).*(n-1) - (7/120)*(Pe^2).*u(n-1)./es ).*Ss./dz;
                 % F42  =  c(n)  - (Pe*Sr + (Pe/es)*Ss)*nts*dt;
                 F4     = [F41; F4; F42];
        end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% linear sink Highest at bottom %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 1
        % U variable
        a1        = ( Pe - (7/45)*(Pe^2).*v(1:n-1) ).*c(1:n-1)./dz ...
                  - (4/105)*((Pe/dz)^2).*u.*(c(2:n) - c(1:n-1)) ...
                  - 2*(7/120)*(Ss/es)*(Pe^2)*dz.*((1:n-1)')/dz;
        a2        = - ( Pe - (7/45)*(Pe^2).*v(1:n-2) ).*c(1:n-2)./dz ...
                  + (4/105)*((Pe/dz)^2).*u(1:n-2).*(c(2:n-1) - c(1:n-2)) ...
                  + 2*(7/120)*(Ss/es)*(Pe^2)*dz.*((1:n-2)')/dz;
        CU        = diag(a1) + diag(a2,-1);
        CU        = [CU; zeros(1,n-1)];
        CU(n,n-1) =  - ( Pe - (7/45)*(Pe^2).*v(n-1) ).*c(n-1)/dz ...
                  + (4/105)*((Pe/dz)^2).*u(n-1).*(c(n) - c(n-1)) ...
                  + 2*(7/120)*(Ss/es)*(Pe^2)*dz*(n-1)/dz;
        % V variable
        a3        = - (7/45)*(Pe^2).*u(1:n-1).*c(1:n-1)./dz;
        a4        = (7/45)*(Pe^2).*u(1:n-2).*c(1:n-2)./dz;
        CV        = diag(a3) + diag(a4,-1);
        CV        = [CV; zeros(1,n-1)];
        CV        = [CV, zeros(n,1)];
        CV(n,n-1) = (7/45)*(Pe^2).*u(n-1).*c(n-1)./dz;
        CV(1,1)   = 0;

        switch s
            case 1
                 % C variable
                 diagC  = 1/dt + ( Pe - (7/45)*(Pe^2).*v(2:n-1) ).*u(2:n-1)./dz ...
                        + ( (2/105)*(Pe^2).*( u(2:n-1).^2 ) + es^2 )./(dz^2)...
                        + ( (2/105)*(Pe^2).*( u(1:n-2).^2 ) + es^2 )./(dz^2);
                 diagC1 = 1/dt + Pe*u(1)/dz ...
                        + ( (2/105)*(Pe^2).*( u(1).^2 ) + es^2 )./(dz^2);
                 diagC2 = 1/dt + ( (2/105)*(Pe^2).*( u(n-1).^2 ) + es^2 )./(dz^2);
                 diagC  = [diagC1; diagC; diagC2];
                 uppC   = - ( (2/105)*(Pe^2).*( u(1:n-1).^2 ) + es^2 )./(dz^2);
                 lowC   = - ( Pe - (7/45)*(Pe^2).*v(1:n-1) ).*u./dz ...
                        - ( (2/105)*(Pe^2).*( u.^2 ) + es^2 )./(dz^2);
                 CC     = diag(diagC) + diag(uppC,1) + diag(lowC,-1);
                 % Nu variable
                 CN     = zeros(n);
                 % Equality
                 F4     = (c(2:n-1) - cp(2:n-1))./dt ...
                        + (Pe - (7/45)*(Pe^2).*v(2:n-1)).*u(2:n-1).*c(2:n-1)./dz...
                        - ( (2/105)*(Pe^2).*(u(2:n-1).^2) + es^2 ).*(c(3:n) - c(2:n-1))./(dz^2)...
                        + 2.*( (Pe*(dz^2)/(2*es)).*((2:n-1)'.^2) ...
                        - (7/120)*(Pe^2)*dz.*u(2:n-1).*((2:n-1)')./es ).*Ss./dz...
                        - (Pe - (7/45)*(Pe^2).*v(1:n-2)).*u(1:n-2).*c(1:n-2)./dz...
                        + ( (2/105)*(Pe^2).*(u(1:n-2).^2) + es^2 ).*(c(2:n-1) - c(1:n-2))./(dz^2)...
                        - 2.*( (Pe*(dz^2)/(2*es)).*((1:n-2)'.^2) ...
                        - (7/120)*(Pe^2)*dz.*u(1:n-2).*((1:n-2)')./es ).*Ss./dz;
                 F41    = (c(1) - cp(1))./dt ...
                        + (Pe - (7/45)*(Pe^2).*v(1)).*u(1).*c(1)./dz...
                        - ( (2/105)*(Pe^2).*(u(1).^2) + es^2 ).*(c(2) - c(1))./(dz^2)...
                        + 2*( (Pe*(dz^2)/(2*es)) - (7/120)*(Pe^2)*dz.*u(1)./es ).*Ss./dz...
                        - Pe*Sl/dz;
                 F42    = (c(n) - cp(n))./dt ...
                        + (Pe*Sr + (Pe/es)*Ss)/dz ...
                        - (Pe - (7/45)*(Pe^2).*v(n-1)).*u(n-1).*c(n-1)./dz...
                        + ( (2/105)*(Pe^2).*(u(n-1).^2) + es^2 ).*(c(n) - c(n-1))./(dz^2)...
                        - 2*( (Pe*(dz^2)/(2*es)).*((n-1)^2) ...
                        - (7/120)*(Pe^2)*dz.*u(n-1)*(n-1)./es ).*Ss./dz;
                 F4     = [F41; F4; F42];


            case 0
                % C variable
                 diagC  = ( Pe - (7/45)*(Pe^2).*v(2:n-1) ).*u(2:n-1)./dz ...
                        + ( (2/105)*(Pe^2).*( u(2:n-1).^2 ) + es^2 )./(dz^2)...
                        + ( (2/105)*(Pe^2).*( u(1:n-2).^2 ) + es^2 )./(dz^2);
                 diagC1 = Pe*u(1)/dz ...
                        + ( (2/105)*(Pe^2).*( u(1).^2 ) + es^2 )./(dz^2);
                 diagC2 = ( (2/105)*(Pe^2).*( u(n-1).^2 ) + es^2 )./(dz^2);
                 % diagC2 = 1;
                 diagC  = [diagC1; diagC; diagC2];
                 uppC   = - ( (2/105)*(Pe^2).*( u(1:n-1).^2 ) + es^2 )./(dz^2);
                 lowC   = - ( Pe - (7/45)*(Pe^2).*v(1:n-1) ).*u./dz ...
                        - ( (2/105)*(Pe^2).*( u.^2 ) + es^2 )./(dz^2);
                 CC     = diag(diagC) + diag(uppC,1) + diag(lowC,-1);
                 % Nu variable
                 CN     = zeros(n);
                 % Equality
                 F4     = (Pe - (7/45)*(Pe^2).*v(2:n-1)).*u(2:n-1).*c(2:n-1)./dz...
                        - ( (2/105)*(Pe^2).*(u(2:n-1).^2) + es^2 ).*(c(3:n) - c(2:n-1))./(dz^2)...
                        + 2.*( (Pe*(dz^2)/(2*es)).*((2:n-1)'.^2) ...
                        - (7/120)*(Pe^2)*dz.*u(2:n-1).*((2:n-1)')./es ).*Ss./dz...
                        - (Pe - (7/45)*(Pe^2).*v(1:n-2)).*u(1:n-2).*c(1:n-2)./dz...
                        + ( (2/105)*(Pe^2).*(u(1:n-2).^2) + es^2 ).*(c(2:n-1) - c(1:n-2))./(dz^2)...
                        - 2.*( (Pe*(dz^2)/(2*es)).*((1:n-2)'.^2) ...
                        - (7/120)*(Pe^2)*dz.*u(1:n-2).*((1:n-2)')./es ).*Ss./dz;
                 F41    = (Pe - (7/45)*(Pe^2).*v(1)).*u(1).*c(1)./dz...
                        - ( (2/105)*(Pe^2).*(u(1).^2) + es^2 ).*(c(2) - c(1))./(dz^2)...
                        + 2*( (Pe*(dz^2)/(2*es)) - (7/120)*(Pe^2)*dz.*u(1)./es ).*Ss./dz...
                        - Pe*Sl/dz;
                 F42    = (Pe*Sr + (Pe/es)*Ss)/dz ...
                        - (Pe - (7/45)*(Pe^2).*v(n-1)).*u(n-1).*c(n-1)./dz...
                        + ( (2/105)*(Pe^2).*(u(n-1).^2) + es^2 ).*(c(n) - c(n-1))./(dz^2)...
                        - 2*( (Pe*(dz^2)/(2*es)).*((n-1)^2) ...
                        - (7/120)*(Pe^2)*dz.*u(n-1)*(n-1)./es ).*Ss./dz;
                 % F42  =  c(n)  - (Pe*Sr + (Pe/es)*Ss)*nts*dt;
                 F4     = [F41; F4; F42];
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% linear sink Highest at top %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 2
        % U variable
        a1        = ( Pe - (7/45)*(Pe^2).*v(1:n-1) ).*c(1:n-1)./dz ...
                  - (4/105)*((Pe/dz)^2).*u.*(c(2:n) - c(1:n-1)) ...
                  - 2*(7/120)*(Ss/es)*(Pe^2).*(1 - dz.*(1:n-1)')/dz;
        a2        = - ( Pe - (7/45)*(Pe^2).*v(1:n-2) ).*c(1:n-2)./dz ...
                  + (4/105)*((Pe/dz)^2).*u(1:n-2).*(c(2:n-1) - c(1:n-2)) ...
                  + 2*(7/120)*(Ss/es)*(Pe^2).*(1 - dz.*(1:n-2)')/dz;
        CU        = diag(a1) + diag(a2,-1);
        CU        = [CU; zeros(1,n-1)];
        CU(n,n-1) =  - ( Pe - (7/45)*(Pe^2).*v(n-1) ).*c(n-1)/dz ...
                  + (4/105)*((Pe/dz)^2).*u(n-1).*(c(n) - c(n-1)) ...
                  + 2*(7/120)*(Ss/es)*(1 - dz*(n-1))*(Pe^2)/dz;
        % V variable
        a3        = - (7/45)*(Pe^2).*u(1:n-1).*c(1:n-1)./dz;
        a4        = (7/45)*(Pe^2).*u(1:n-2).*c(1:n-2)./dz;
        CV        = diag(a3) + diag(a4,-1);
        CV        = [CV; zeros(1,n-1)];
        CV        = [CV, zeros(n,1)];
        CV(n,n-1) = (7/45)*(Pe^2).*u(n-1).*c(n-1)./dz;
        CV(1,1)   = 0;

        switch s
            case 1
                 % C variable
                 diagC  = 1/dt + ( Pe - (7/45)*(Pe^2).*v(2:n-1) ).*u(2:n-1)./dz ...
                        + ( (2/105)*(Pe^2).*( u(2:n-1).^2 ) + es^2 )./(dz^2)...
                        + ( (2/105)*(Pe^2).*( u(1:n-2).^2 ) + es^2 )./(dz^2);
                 diagC1 = 1/dt + Pe*u(1)/dz ...
                        + ( (2/105)*(Pe^2).*( u(1).^2 ) + es^2 )./(dz^2);
                 diagC2 = 1/dt + ( (2/105)*(Pe^2).*( u(n-1).^2 ) + es^2 )./(dz^2);
                 diagC  = [diagC1; diagC; diagC2];
                 uppC   = - ( (2/105)*(Pe^2).*( u(1:n-1).^2 ) + es^2 )./(dz^2);
                 lowC   = - ( Pe - (7/45)*(Pe^2).*v(1:n-1) ).*u./dz ...
                        - ( (2/105)*(Pe^2).*( u.^2 ) + es^2 )./(dz^2);
                 CC     = diag(diagC) + diag(uppC,1) + diag(lowC,-1);
                 % Nu variable
                 CN     = zeros(n);
                 % Equality
                 F4     = (c(2:n-1) - cp(2:n-1))./dt ...
                        + (Pe - (7/45)*(Pe^2).*v(2:n-1)).*u(2:n-1).*c(2:n-1)./dz...
                        - ( (2/105)*(Pe^2).*(u(2:n-1).^2) + es^2 ).*(c(3:n) - c(2:n-1))./(dz^2)...
                        + 2.*( (Pe*dz/es).*((2:n-1)').*(1 - dz.*((2:n-1)')./2) ...
                        - (7/120)*(Pe^2).*u(2:n-1).*(1 - dz.*(2:n-1)')./es ).*Ss./dz...
                        - (Pe - (7/45)*(Pe^2).*v(1:n-2)).*u(1:n-2).*c(1:n-2)./dz...
                        + ( (2/105)*(Pe^2).*(u(1:n-2).^2) + es^2 ).*(c(2:n-1) - c(1:n-2))./(dz^2)...
                        - 2.*( (Pe*dz/es).*((1:n-2)').*(1 - dz.*(1:n-2)'./2) ...
                        - (7/120)*(Pe^2).*u(1:n-2).*(1 - dz.*(1:n-2)')./es ).*Ss./dz;
                 F41    = (c(1) - cp(1))./dt ...
                        + (Pe - (7/45)*(Pe^2).*v(1)).*u(1).*c(1)./dz...
                        - ( (2/105)*(Pe^2).*(u(1).^2) + es^2 ).*(c(2) - c(1))./(dz^2)...
                        + 2*( (Pe*dz/es)*(1 - dz/2) - (7/120)*(Pe^2).*u(1)*(1 - dz)./es ).*Ss./dz...
                        - Pe*Sl/dz;
                 F42    = (c(n) - cp(n))./dt ...
                        + (Pe*Sr + (Pe/es)*Ss)/dz ...
                        - (Pe - (7/45)*(Pe^2).*v(n-1)).*u(n-1).*c(n-1)./dz...
                        + ( (2/105)*(Pe^2).*(u(n-1).^2) + es^2 ).*(c(n) - c(n-1))./(dz^2)...
                        - 2*( (Pe*dz/es)*(n-1)*(1 - dz*(n-1)/2) ...
                        - (7/120)*(Pe^2).*u(n-1)*(1 - dz*(n-1))./es ).*Ss./dz;
                 F4     = [F41; F4; F42];

            case 0
                % C variable
                 diagC  = ( Pe - (7/45)*(Pe^2).*v(2:n-1) ).*u(2:n-1)./dz ...
                        + ( (2/105)*(Pe^2).*( u(2:n-1).^2 ) + es^2 )./(dz^2)...
                        + ( (2/105)*(Pe^2).*( u(1:n-2).^2 ) + es^2 )./(dz^2);
                 diagC1 = Pe*u(1)/dz ...
                        + ( (2/105)*(Pe^2).*( u(1).^2 ) + es^2 )./(dz^2);
                 diagC2 = ( (2/105)*(Pe^2).*( u(n-1).^2 ) + es^2 )./(dz^2);
                 % diagC2 = 1;
                 diagC  = [diagC1; diagC; diagC2];
                 uppC   = - ( (2/105)*(Pe^2).*( u(1:n-1).^2 ) + es^2 )./(dz^2);
                 lowC   = - ( Pe - (7/45)*(Pe^2).*v(1:n-1) ).*u./dz ...
                        - ( (2/105)*(Pe^2).*( u.^2 ) + es^2 )./(dz^2);
                 CC     = diag(diagC) + diag(uppC,1) + diag(lowC,-1);
                 % Nu variable
                 CN     = zeros(n);
                 % Equality
                 F4     = (Pe - (7/45)*(Pe^2).*v(2:n-1)).*u(2:n-1).*c(2:n-1)./dz...
                        - ( (2/105)*(Pe^2).*(u(2:n-1).^2) + es^2 ).*(c(3:n) - c(2:n-1))./(dz^2)...
                        + 2.*( (Pe*dz/es).*((2:n-1)').*(1 - dz.*((2:n-1)')./2) ...
                        - (7/120)*(Pe^2).*u(2:n-1).*(1 - dz.*(2:n-1)')./es ).*Ss./dz...
                        - (Pe - (7/45)*(Pe^2).*v(1:n-2)).*u(1:n-2).*c(1:n-2)./dz...
                        + ( (2/105)*(Pe^2).*(u(1:n-2).^2) + es^2 ).*(c(2:n-1) - c(1:n-2))./(dz^2)...
                        - 2.*( (Pe*dz/es).*((1:n-2)').*(1 - dz.*(1:n-2)'./2) ...
                        - (7/120)*(Pe^2).*u(1:n-2).*(1 - dz.*(1:n-2)')./es ).*Ss./dz;
                 F41    = (Pe - (7/45)*(Pe^2).*v(1)).*u(1).*c(1)./dz...
                        - ( (2/105)*(Pe^2).*(u(1).^2) + es^2 ).*(c(2) - c(1))./(dz^2)...
                        + 2*( (Pe*dz/es)*(1 - dz/2) - (7/120)*(Pe^2).*u(1)*(1 - dz)./es ).*Ss./dz...
                        - Pe*Sl/dz;
                 F42    = (Pe*Sr + (Pe/es)*Ss)/dz ...
                        - (Pe - (7/45)*(Pe^2).*v(n-1)).*u(n-1).*c(n-1)./dz...
                        + ( (2/105)*(Pe^2).*(u(n-1).^2) + es^2 ).*(c(n) - c(n-1))./(dz^2)...
                        - 2*( (Pe*dz/es)*(n-1)*(1 - dz*(n-1)/2) ...
                        - (7/120)*(Pe^2).*u(n-1)*(1 - dz*(n-1))./es ).*Ss./dz;
                 % F42  =  c(n)  - (Pe*Sr + (Pe/es)*Ss)*nts*dt;
                 F4     = [F41; F4; F42];
        end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2nd order polynomial sink Highest at bottom %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 3
        % U variable
        a1        = ( Pe - (7/45)*(Pe^2).*v(1:n-1) ).*c(1:n-1)./dz ...
                  - (4/105)*((Pe/dz)^2).*u.*(c(2:n) - c(1:n-1)) ...
                  - 3*(7/120)*(Ss/es)*(Pe^2)*(dz^2).*((1:n-1)'.^2)/dz;
        a2        = - ( Pe - (7/45)*(Pe^2).*v(1:n-2) ).*c(1:n-2)./dz ...
                  + (4/105)*((Pe/dz)^2).*u(1:n-2).*(c(2:n-1) - c(1:n-2)) ...
                  + 3*(7/120)*(Ss/es)*(Pe^2)*(dz^2).*((1:n-2)'.^2)/dz;
        CU        = diag(a1) + diag(a2,-1);
        CU        = [CU; zeros(1,n-1)];
        CU(n,n-1) =  - ( Pe - (7/45)*(Pe^2).*v(n-1) ).*c(n-1)/dz ...
                  + (4/105)*((Pe/dz)^2).*u(n-1).*(c(n) - c(n-1)) ...
                  + 3*(7/120)*(Ss/es)*(dz^2)*((n-1)^2)*(Pe^2)/dz;
        % V variable
        a3        = - (7/45)*(Pe^2).*u(1:n-1).*c(1:n-1)./dz;
        a4        = (7/45)*(Pe^2).*u(1:n-2).*c(1:n-2)./dz;
        CV        = diag(a3) + diag(a4,-1);
        CV        = [CV; zeros(1,n-1)];
        CV        = [CV, zeros(n,1)];
        CV(n,n-1) = (7/45)*(Pe^2).*u(n-1).*c(n-1)./dz;
        CV(1,1)   = 0;

        switch s
            case 1
                 % C variable
                 diagC  = 1/dt + ( Pe - (7/45)*(Pe^2).*v(2:n-1) ).*u(2:n-1)./dz ...
                        + ( (2/105)*(Pe^2).*( u(2:n-1).^2 ) + es^2 )./(dz^2)...
                        + ( (2/105)*(Pe^2).*( u(1:n-2).^2 ) + es^2 )./(dz^2);
                 diagC1 = 1/dt + Pe*u(1)/dz ...
                        + ( (2/105)*(Pe^2).*( u(1).^2 ) + es^2 )./(dz^2);
                 diagC2 = 1/dt + ( (2/105)*(Pe^2).*( u(n-1).^2 ) + es^2 )./(dz^2);
                 diagC  = [diagC1; diagC; diagC2];
                 uppC   = - ( (2/105)*(Pe^2).*( u(1:n-1).^2 ) + es^2 )./(dz^2);
                 lowC   = - ( Pe - (7/45)*(Pe^2).*v(1:n-1) ).*u./dz ...
                        - ( (2/105)*(Pe^2).*( u.^2 ) + es^2 )./(dz^2);
                 CC     = diag(diagC) + diag(uppC,1) + diag(lowC,-1);
                 % Nu variable
                 CN     = zeros(n);
                 % Equality
                 F4     = (c(2:n-1) - cp(2:n-1))./dt ...
                        + (Pe - (7/45)*(Pe^2).*v(2:n-1)).*u(2:n-1).*c(2:n-1)./dz...
                        - ( (2/105)*(Pe^2).*(u(2:n-1).^2) + es^2 ).*(c(3:n) - c(2:n-1))./(dz^2)...
                        + 3.*( (Pe*(dz^3)/es).*((2:n-1)'.^3)./3 ...
                        - (7/120)*(Pe^2)*(dz^2).*u(2:n-1).*((2:n-1)'.^2)./es ).*Ss./dz...
                        - (Pe - (7/45)*(Pe^2).*v(1:n-2)).*u(1:n-2).*c(1:n-2)./dz...
                        + ( (2/105)*(Pe^2).*(u(1:n-2).^2) + es^2 ).*(c(2:n-1) - c(1:n-2))./(dz^2)...
                        - 3.*( (Pe*(dz^3)/es).*((1:n-2)'.^3)./3 ...
                        - (7/120)*(Pe^2)*(dz^2).*u(1:n-2).*((1:n-2)'.^2)./es ).*Ss./dz;
                 F41    = (c(1) - cp(1))./dt ...
                        + (Pe - (7/45)*(Pe^2).*v(1)).*u(1).*c(1)./dz...
                        - ( (2/105)*(Pe^2).*(u(1).^2) + es^2 ).*(c(2) - c(1))./(dz^2)...
                        + 3*( (Pe*(dz^3)/es)/3 - (7/120)*(dz^2)*(Pe^2)*u(1)/es ).*Ss./dz...
                        - Pe*Sl/dz;
                 F42    = (c(n) - cp(n))./dt ...
                        + (Pe*Sr + (Pe/es)*Ss)/dz ...
                        - (Pe - (7/45)*(Pe^2).*v(n-1)).*u(n-1).*c(n-1)./dz...
                        + ( (2/105)*(Pe^2).*(u(n-1).^2) + es^2 ).*(c(n) - c(n-1))./(dz^2)...
                        - 3*( (Pe*(dz^3)/es)*(n-1)/3 ...
                        - (7/120)*(Pe^2).*u(n-1)*(dz^2)*((n-1)^2)./es ).*Ss./dz;
                 F4     = [F41; F4; F42];

            case 0
                % C variable
                 diagC  = ( Pe - (7/45)*(Pe^2).*v(2:n-1) ).*u(2:n-1)./dz ...
                        + ( (2/105)*(Pe^2).*( u(2:n-1).^2 ) + es^2 )./(dz^2)...
                        + ( (2/105)*(Pe^2).*( u(1:n-2).^2 ) + es^2 )./(dz^2);
                 diagC1 = Pe*u(1)/dz ...
                        + ( (2/105)*(Pe^2).*( u(1).^2 ) + es^2 )./(dz^2);
                 diagC2 = ( (2/105)*(Pe^2).*( u(n-1).^2 ) + es^2 )./(dz^2);
                 % diagC2 = 1;
                 diagC  = [diagC1; diagC; diagC2];
                 uppC   = - ( (2/105)*(Pe^2).*( u(1:n-1).^2 ) + es^2 )./(dz^2);
                 lowC   = - ( Pe - (7/45)*(Pe^2).*v(1:n-1) ).*u./dz ...
                        - ( (2/105)*(Pe^2).*( u.^2 ) + es^2 )./(dz^2);
                 CC     = diag(diagC) + diag(uppC,1) + diag(lowC,-1);
                 % Nu variable
                 CN     = zeros(n);
                 % Equality
                 F4     = (Pe - (7/45)*(Pe^2).*v(2:n-1)).*u(2:n-1).*c(2:n-1)./dz...
                        - ( (2/105)*(Pe^2).*(u(2:n-1).^2) + es^2 ).*(c(3:n) - c(2:n-1))./(dz^2)...
                        + 3.*( (Pe*(dz^3)/es).*((2:n-1)'.^3)./3 ...
                        - (7/120)*(Pe^2)*(dz^2).*u(2:n-1).*((2:n-1)'.^2)./es ).*Ss./dz...
                        - (Pe - (7/45)*(Pe^2).*v(1:n-2)).*u(1:n-2).*c(1:n-2)./dz...
                        + ( (2/105)*(Pe^2).*(u(1:n-2).^2) + es^2 ).*(c(2:n-1) - c(1:n-2))./(dz^2)...
                        - 3.*( (Pe*(dz^3)/es).*((1:n-2)'.^3)./3 ...
                        - (7/120)*(Pe^2)*(dz^2).*u(1:n-2).*((1:n-2)'.^2)./es ).*Ss./dz;
                 F41    = (Pe - (7/45)*(Pe^2).*v(1)).*u(1).*c(1)./dz...
                        - ( (2/105)*(Pe^2).*(u(1).^2) + es^2 ).*(c(2) - c(1))./(dz^2)...
                        + 3*( (Pe*(dz^3)/es)/3 - (7/120)*(dz^2)*(Pe^2)*u(1)/es ).*Ss./dz...
                        - Pe*Sl/dz;
                 F42    = (Pe*Sr + (Pe/es)*Ss)/dz ...
                        - (Pe - (7/45)*(Pe^2).*v(n-1)).*u(n-1).*c(n-1)./dz...
                        + ( (2/105)*(Pe^2).*(u(n-1).^2) + es^2 ).*(c(n) - c(n-1))./(dz^2)...
                        - 3*( (Pe*(dz^3)/es)*(n-1)/3 ...
                        - (7/120)*(Pe^2).*u(n-1)*(dz^2)*((n-1)^2)./es ).*Ss./dz;
                 % F42  =  c(n)  - (Pe*Sr + (Pe/es)*Ss)*nts*dt;
                 F4     = [F41; F4; F42];
        end


%%%%%%%%%%%%%%%%%%%%%%%%% cts respiration and linear growth highest at top %%%%%%%%%%%%%%%%%%%%%%%%
    case 4
        % U variable
        a1        = ( Pe - (7/45)*(Pe^2).*v(1:n-1) ).*c(1:n-1)./dz ...
                  - (4/105)*((Pe/dz)^2).*u.*(c(2:n) - c(1:n-1)) ...
                  - 2.*(7/120)*(Ss2/es)*(Pe^2).*( dz.*(1:n-1)')./dz ...
                  + (7/120)*(Ss1/es)*(Pe^2)/dz;
        a2        = - ( Pe - (7/45)*(Pe^2).*v(1:n-2) ).*c(1:n-2)./dz ...
                  + (4/105)*((Pe/dz)^2).*u(1:n-2).*(c(2:n-1) - c(1:n-2)) ...
                  + 2.*(7/120)*(Ss2/es)*(Pe^2).*( dz.*(1:n-2)')./dz ...
                  + (7/120)*(Ss1/es)*(Pe^2)/dz;
        CU        = diag(a1) + diag(a2,-1);
        CU        = [CU; zeros(1,n-1)];
        CU(n,n-1) =  - ( Pe - (7/45)*(Pe^2).*v(n-1) ).*c(n-1)/dz ...
                  + (4/105)*((Pe/dz)^2).*u(n-1).*(c(n) - c(n-1)) ...
                  + 2.*(7/120)*(Ss2/es)*(Pe^2).*(dz.*(n-1))/dz ...
                  + (7/120)*(Ss1/es)*(Pe^2)/dz;
        % V variable
        a3        = - (7/45)*(Pe^2).*u(1:n-1).*c(1:n-1)./dz;
        a4        = (7/45)*(Pe^2).*u(1:n-2).*c(1:n-2)./dz;
        CV        = diag(a3) + diag(a4,-1);
        CV        = [CV; zeros(1,n-1)];
        CV        = [CV, zeros(n,1)];
        CV(n,n-1) = (7/45)*(Pe^2).*u(n-1).*c(n-1)./dz;
        CV(1,1)   = 0;

        switch s
            case 1
                 % C variable
                 diagC  = 1/dt + ( Pe - (7/45)*(Pe^2).*v(2:n-1) ).*u(2:n-1)./dz ...
                        + ( (2/105)*(Pe^2).*( u(2:n-1).^2 ) + es^2 )./(dz^2)...
                        + ( (2/105)*(Pe^2).*( u(1:n-2).^2 ) + es^2 )./(dz^2);
                 diagC1 = 1/dt + Pe*u(1)/dz ...
                        + ( (2/105)*(Pe^2).*( u(1).^2 ) + es^2 )./(dz^2);
                 diagC2 = 1/dt + ( (2/105)*(Pe^2).*( u(n-1).^2 ) + es^2 )./(dz^2);
                 diagC  = [diagC1; diagC; diagC2];
                 uppC   = - ( (2/105)*(Pe^2).*( u(1:n-1).^2 ) + es^2 )./(dz^2);
                 lowC   = - ( Pe - (7/45)*(Pe^2).*v(1:n-1) ).*u./dz ...
                        - ( (2/105)*(Pe^2).*( u.^2 ) + es^2 )./(dz^2);
                 CC     = diag(diagC) + diag(uppC,1) + diag(lowC,-1);
                 % Nu variable
                 CN     = zeros(n);
                 % Equality
                 F4     = (c(2:n-1) - cp(2:n-1))./dt ...
                        + (Pe - (7/45)*(Pe^2).*v(2:n-1)).*u(2:n-1).*c(2:n-1)./dz...
                        - ( (2/105)*(Pe^2).*(u(2:n-1).^2) + es^2 ).*(c(3:n) - c(2:n-1))./(dz^2)...
                        + ( (Pe*dz/es).*((2:n-1)') - (7/120)*( (Pe^2)/es ).*u(2:n-1) ).*Ss1./dz ...
                        + ( (Pe*(dz^2)/es).*((2:n-1)'.^2) - (7/60)*( (Pe^2)/es )*dz.*u(2:n-1).*((2:n-1)') ).*Ss2./dz ...
                        - (Pe - (7/45)*(Pe^2).*v(1:n-2)).*u(1:n-2).*c(1:n-2)./dz...
                        + ( (2/105)*(Pe^2).*(u(1:n-2).^2) + es^2 ).*(c(2:n-1) - c(1:n-2))./(dz^2)...
                        - ( (Pe*dz/es).*((1:n-2)') - (7/120)*( (Pe^2)/es ).*u(1:n-2) ).*Ss1./dz ...
                        - ( (Pe*(dz^2)/es).*((1:n-2)'.^2) - (7/60)*( (Pe^2)/es )*dz.*u(1:n-2).*((1:n-2)') ).*Ss2./dz;
                 F41    = (c(1) - cp(1))./dt ...
                        + (Pe - (7/45)*(Pe^2).*v(1)).*u(1).*c(1)./dz...
                        - ( (2/105)*(Pe^2).*(u(1).^2) + es^2 ).*(c(2) - c(1))./(dz^2)...
                        + ( (Pe*dz/es) - (7/120)*( (Pe^2)/es ).*u(1) ).*Ss1/dz ...
                        + ( (Pe*(dz^2)/es) - (7/60)*( (Pe^2)/es )*dz.*u(1) )*Ss2/dz ...
                        - Pe*Sl/dz;
                 F42    = (c(n) - cp(n))./dt ...
                        + (Pe*Sr + (Pe/es)*Ss)/dz ...
                        - (Pe - (7/45)*(Pe^2).*v(n-1)).*u(n-1).*c(n-1)./dz...
                        + ( (2/105)*(Pe^2).*(u(n-1).^2) + es^2 ).*(c(n) - c(n-1))./(dz^2)...
                        - ( (Pe*dz/es).*((n-1)) - (7/120)*( (Pe^2)/es ).*u(n-1) ).*Ss1/dz ...
                        - ( (Pe*(dz^2)/es).*((n-1).^2) - (7/60)*( (Pe^2)/es )*dz.*u(n-1).*((n-1)) ).*Ss2/dz;
                 F4     = [F41; F4; F42];

            case 0
                % C variable
                 diagC  = ( Pe - (7/45)*(Pe^2).*v(2:n-1) ).*u(2:n-1)./dz ...
                        + ( (2/105)*(Pe^2).*( u(2:n-1).^2 ) + es^2 )./(dz^2)...
                        + ( (2/105)*(Pe^2).*( u(1:n-2).^2 ) + es^2 )./(dz^2);
                 diagC1 = Pe*u(1)/dz ...
                        + ( (2/105)*(Pe^2).*( u(1).^2 ) + es^2 )./(dz^2);
                 diagC2 = ( (2/105)*(Pe^2).*( u(n-1).^2 ) + es^2 )./(dz^2);
                 % diagC2 = 1;
                 diagC  = [diagC1; diagC; diagC2];
                 uppC   = - ( (2/105)*(Pe^2).*( u(1:n-1).^2 ) + es^2 )./(dz^2);
                 lowC   = - ( Pe - (7/45)*(Pe^2).*v(1:n-1) ).*u./dz ...
                        - ( (2/105)*(Pe^2).*( u.^2 ) + es^2 )./(dz^2);
                 CC     = diag(diagC) + diag(uppC,1) + diag(lowC,-1);
                 % Nu variable
                 CN     = zeros(n);
                 % Equality
                 F4     = (Pe - (7/45)*(Pe^2).*v(2:n-1)).*u(2:n-1).*c(2:n-1)./dz...
                        - ( (2/105)*(Pe^2).*(u(2:n-1).^2) + es^2 ).*(c(3:n) - c(2:n-1))./(dz^2)...
                        + ( (Pe*dz/es).*((2:n-1)') - (7/120)*( (Pe^2)/es ).*u(2:n-1) ).*Ss1/dz ...
                        + ( (Pe*(dz^2)/es).*((2:n-1)'.^2) - (7/60)*( (Pe^2)/es )*dz.*u(2:n-1).*((2:n-1)') ).*Ss2/dz ...
                        - (Pe - (7/45)*(Pe^2).*v(1:n-2)).*u(1:n-2).*c(1:n-2)./dz...
                        + ( (2/105)*(Pe^2).*(u(1:n-2).^2) + es^2 ).*(c(2:n-1) - c(1:n-2))./(dz^2)...
                        - ( (Pe*dz/es).*((1:n-2)') - (7/120)*( (Pe^2)/es ).*u(1:n-2) ).*Ss1/dz ...
                        - ( (Pe*(dz^2)/es).*((1:n-2)'.^2) - (7/60)*( (Pe^2)/es )*dz.*u(1:n-2).*((1:n-2)') ).*Ss2/dz;
                 F41    = (Pe - (7/45)*(Pe^2).*v(1)).*u(1).*c(1)./dz...
                        - ( (2/105)*(Pe^2).*(u(1).^2) + es^2 ).*(c(2) - c(1))./(dz^2)...
                        + ( (Pe*dz/es) - (7/120)*( (Pe^2)/es ).*u(1) ).*Ss1/dz ...
                        + ( (Pe*(dz^2)/es) - (7/60)*( (Pe^2)/es )*dz.*u(1) ).*Ss2/dz ...
                        - Pe*Sl/dz;
                 F42    = (Pe*Sr + (Pe/es)*Ss)/dz ...
                        - (Pe - (7/45)*(Pe^2)*v(n-1))*u(n-1)*c(n-1)/dz...
                        + ( (2/105)*(Pe^2)*(u(n-1)^2) + es^2 )*(c(n) - c(n-1))/(dz^2)...
                        - ( (Pe*dz/es).*((n-1)) - (7/120)*( (Pe^2)/es )*u(n-1) ).*Ss1/dz ...
                        - ( (Pe*(dz^2)/es).*((n-1)^2) - (7/60)*( (Pe^2)/es )*dz*u(n-1)*((n-1)) ).*Ss2/dz;
                 % F42  =  c(n)  - (Pe*Sr + (Pe/es)*Ss)*nts*dt;
                 F4     = [F41; F4; F42];
        end
              
%%%%%%%%%%%%%%%%%%%%%%%%% cts respiration and linear growth highest at botttom %%%%%%%%%%%%%%%%%%%%%%%%
    case 5
        % U variable
        a1        = ( Pe - (7/45)*(Pe^2).*v(1:n-1) ).*c(1:n-1)./dz ...
                  - (4/105)*((Pe/dz)^2).*u.*(c(2:n) - c(1:n-1)) ...
                  - 2.*(7/120)*(Ss2/es)*(Pe^2).*( 1 - dz.*(1:n-1)' )/dz ...
                  + (7/120)*(Ss1/es)*(Pe^2)/dz;
        a2        = - ( Pe - (7/45)*(Pe^2).*v(1:n-2) ).*c(1:n-2)./dz ...
                  + (4/105)*((Pe/dz)^2).*u(1:n-2).*(c(2:n-1) - c(1:n-2)) ...
                  + 2.*(7/120)*(Ss2/es)*(Pe^2).*(1 - dz.*(1:n-2)')/dz ...
                  + (7/120)*(Ss1/es)*(Pe^2)/dz;
        CU        = diag(a1) + diag(a2,-1);
        CU        = [CU; zeros(1,n-1)];
        CU(n,n-1) =  - ( Pe - (7/45)*(Pe^2).*v(n-1) ).*c(n-1)/dz ...
                  + (4/105)*((Pe/dz)^2).*u(n-1).*(c(n) - c(n-1)) ...
                  + 2.*(7/120)*(Ss2/es)*(Pe^2).*(1 - dz.*(n-1))/dz ...
                  + (7/120)*(Ss1/es)*(Pe^2)/dz;
        % V variable
        a3        = - (7/45)*(Pe^2).*u(1:n-1).*c(1:n-1)./dz;
        a4        = (7/45)*(Pe^2).*u(1:n-2).*c(1:n-2)./dz;
        CV        = diag(a3) + diag(a4,-1);
        CV        = [CV; zeros(1,n-1)];
        CV        = [CV, zeros(n,1)];
        CV(n,n-1) = (7/45)*(Pe^2).*u(n-1).*c(n-1)./dz;
        CV(1,1)   = 0;

        switch s
            case 1
                 % C variable
                 diagC  = 1/dt + ( Pe - (7/45)*(Pe^2).*v(2:n-1) ).*u(2:n-1)./dz ...
                        + ( (2/105)*(Pe^2).*( u(2:n-1).^2 ) + es^2 )./(dz^2)...
                        + ( (2/105)*(Pe^2).*( u(1:n-2).^2 ) + es^2 )./(dz^2);
                 diagC1 = 1/dt + Pe*u(1)/dz ...
                        + ( (2/105)*(Pe^2).*( u(1).^2 ) + es^2 )./(dz^2);
                 diagC2 = 1/dt + ( (2/105)*(Pe^2).*( u(n-1).^2 ) + es^2 )./(dz^2);
                 diagC  = [diagC1; diagC; diagC2];
                 uppC   = - ( (2/105)*(Pe^2).*( u(1:n-1).^2 ) + es^2 )./(dz^2);
                 lowC   = - ( Pe - (7/45)*(Pe^2).*v(1:n-1) ).*u./dz ...
                        - ( (2/105)*(Pe^2).*( u.^2 ) + es^2 )./(dz^2);
                 CC     = diag(diagC) + diag(uppC,1) + diag(lowC,-1);
                 % Nu variable
                 CN     = zeros(n);
                 % Equality
                 F4     = (c(2:n-1) - cp(2:n-1))./dt ...
                        + (Pe - (7/45)*(Pe^2).*v(2:n-1)).*u(2:n-1).*c(2:n-1)./dz...
                        - ( (2/105)*(Pe^2).*(u(2:n-1).^2) + es^2 ).*(c(3:n) - c(2:n-1))./(dz^2)...
                        + ( (Pe*dz/es).*((2:n-1)') - (7/120)*( (Pe^2)/es ).*u(2:n-1) ).*Ss1 ...
                        + 2.*( (Pe*dz/es).*((2:n-1)').*(1 - dz.*((2:n-1)')./2) ...
                        - (7/120)*(Pe^2).*u(2:n-1).*(1 - dz.*(2:n-1)')./es ).*Ss2./dz ...
                        - (Pe - (7/45)*(Pe^2).*v(1:n-2)).*u(1:n-2).*c(1:n-2)./dz...
                        + ( (2/105)*(Pe^2).*(u(1:n-2).^2) + es^2 ).*(c(2:n-1) - c(1:n-2))./(dz^2)...
                        - ( (Pe*dz/es).*((1:n-2)') - (7/120)*( (Pe^2)/es ).*u(1:n-2) ).*Ss1 ...
                        - 2.*( (Pe*dz/es).*((1:n-2)').*(1 - dz.*(1:n-2)'./2) ...
                        - (7/120)*(Pe^2).*u(1:n-2).*(1 - dz.*(1:n-2)')./es ).*Ss2./dz;

                 F41    = (c(1) - cp(1))./dt ...
                        + (Pe - (7/45)*(Pe^2).*v(1)).*u(1).*c(1)./dz...
                        - ( (2/105)*(Pe^2).*(u(1).^2) + es^2 ).*(c(2) - c(1))./(dz^2)...
                        + ( (Pe*dz/es) - (7/120)*( (Pe^2)/es ).*u(1) ).*Ss1 ...
                        + 2.*( (Pe*dz/es)*(1 - dz/2) - (7/120)*(Pe^2).*u(1).*(1 - dz)./es ).*Ss2./dz ...
                        - Pe*Sl/dz;
                 F42    = (c(n) - cp(n))./dt ...
                        + (Pe*Sr + (Pe/es)*Ss)/dz ...
                        - (Pe - (7/45)*(Pe^2).*v(n-1)).*u(n-1).*c(n-1)./dz...
                        + ( (2/105)*(Pe^2).*(u(n-1).^2) + es^2 ).*(c(n) - c(n-1))./(dz^2)...
                        - ( (Pe*dz/es).*((n-1)) - (7/120)*( (Pe^2)/es ).*u(n-1) ).*Ss1 ...
                        - 2.*( (Pe*dz/es).*((n-1)).*(1 - dz*(n-1)/2) ...
                        - (7/120)*(Pe^2).*u(n-1).*(1 - dz.*(n-1))./es ).*Ss2./dz;
                 F4     = [F41; F4; F42];

            case 0
                % C variable
                 diagC  = ( Pe - (7/45)*(Pe^2).*v(2:n-1) ).*u(2:n-1)./dz ...
                        + ( (2/105)*(Pe^2).*( u(2:n-1).^2 ) + es^2 )./(dz^2)...
                        + ( (2/105)*(Pe^2).*( u(1:n-2).^2 ) + es^2 )./(dz^2);
                 diagC1 = Pe*u(1)/dz ...
                        + ( (2/105)*(Pe^2).*( u(1).^2 ) + es^2 )./(dz^2);
                 diagC2 = ( (2/105)*(Pe^2).*( u(n-1).^2 ) + es^2 )./(dz^2);
                 % diagC2 = 1;
                 diagC  = [diagC1; diagC; diagC2];
                 uppC   = - ( (2/105)*(Pe^2).*( u(1:n-1).^2 ) + es^2 )./(dz^2);
                 lowC   = - ( Pe - (7/45)*(Pe^2).*v(1:n-1) ).*u./dz ...
                        - ( (2/105)*(Pe^2).*( u.^2 ) + es^2 )./(dz^2);
                 CC     = diag(diagC) + diag(uppC,1) + diag(lowC,-1);
                 % Nu variable
                 CN     = zeros(n);
                 % Equality
                 F4     = (Pe - (7/45)*(Pe^2).*v(2:n-1)).*u(2:n-1).*c(2:n-1)./dz...
                        - ( (2/105)*(Pe^2).*(u(2:n-1).^2) + es^2 ).*(c(3:n) - c(2:n-1))./(dz^2)...
                        + ( (Pe*dz/es).*((2:n-1)') - (7/120)*( (Pe^2)/es ).*u(2:n-1) ).*Ss1 ...
                        + 2.*( (Pe*dz/es).*((2:n-1)').*(1 - dz.*((2:n-1)')./2) ...
                        - (7/120)*(Pe^2).*u(2:n-1).*(1 - dz.*(2:n-1)')./es ).*Ss2./dz ...
                        - (Pe - (7/45)*(Pe^2).*v(1:n-2)).*u(1:n-2).*c(1:n-2)./dz...
                        + ( (2/105)*(Pe^2).*(u(1:n-2).^2) + es^2 ).*(c(2:n-1) - c(1:n-2))./(dz^2)...
                        - ( (Pe*dz/es).*((1:n-2)') - (7/120)*( (Pe^2)/es ).*u(1:n-2) ).*Ss1 ...
                        - 2.*( (Pe*dz/es).*((1:n-2)').*(1 - dz.*(1:n-2)'./2) ...
                        - (7/120)*(Pe^2).*u(1:n-2).*(1 - dz.*(1:n-2)')./es ).*Ss2./dz;

                 F41    = (Pe - (7/45)*(Pe^2).*v(1)).*u(1).*c(1)./dz...
                        - ( (2/105)*(Pe^2).*(u(1).^2) + es^2 ).*(c(2) - c(1))./(dz^2)...
                        + ( (Pe*dz/es) - (7/120)*( (Pe^2)/es ).*u(1) ).*Ss1 ...
                        + 2.*( (Pe*dz/es)*(1 - dz/2) - (7/120)*(Pe^2).*u(1).*(1 - dz)./es ).*Ss2./dz ...
                        - Pe*Sl/dz;
                 F42    = (Pe*Sr + (Pe/es)*Ss)/dz ...
                        - (Pe - (7/45)*(Pe^2).*v(n-1)).*u(n-1).*c(n-1)./dz...
                        + ( (2/105)*(Pe^2).*(u(n-1).^2) + es^2 ).*(c(n) - c(n-1))./(dz^2)...
                        - ( (Pe*dz/es).*((n-1)) - (7/120)*( (Pe^2)/es ).*u(n-1) ).*Ss1 ...
                        - 2.*( (Pe*dz/es).*((n-1)).*(1 - dz*(n-1)/2) ...
                        - (7/120)*(Pe^2).*u(n-1).*(1 - dz.*(n-1))./es ).*Ss2./dz;
                 % F42  =  c(n)  - (Pe*Sr + (Pe/es)*Ss)*nts*dt;
                 F4     = [F41; F4; F42];
        end

end



end
