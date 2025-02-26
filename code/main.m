% clear all
close all
addpath('./fcn/');
addpath('./solver/');
addpath('./initial/');
path = "/Users/mazennakad/Desktop/phloem/phloem_biomeE";
data = readtable(path + "/output/BCI_hydro_Cohort_hourly.csv");
% duration = height(data)/600;
duration = 24;

% Get a c-nu relation at ambient temperature
c = linspace(0,2000,1000)';
[~,bb,rrs] = MLR_Viscosity(c,293.15,1000);

%% Grid and time step
n   = 100;   % nb of cells      (at least 300 for accuracy)
dz  = 1/n;  % size of the cell
dt  = 1e1;  % timestep
nts = 1440; % nb of time steps (this makes it 1 hour)

%% Constants
k     = 5e-14;                       % membrane permeability m/(Pa s)
a     = 1e-5;                        % phloem thickness m
T     = 293.15;                      % temperature K
Rg    = 8.3145;                      % Ideal gas constant J/(mol K)
Mw    = 342.2965;                    % Molecular mass of sucrose g/mol
rho   = 1e3;                         % Density of water kg/m3
g     = 9.81;                        % Gravitational constant m/s2
D     = 4e-10;                       % Sucrose diffusivity in water m2/s

%%% Declare variables
c    = zeros(n,nts);   % Sucrose concentration at the center of the cell
u    = zeros(n-1,nts); % Axial velocity at the face of the cell
p    = zeros(n,nts);   % Dynamic fluid pressure at the center of the cell
v    = zeros(n,nts);   % Radial velocity at the center of the cell
nu   = zeros(n,nts);   % viscosity at the center of the cell

% cs1         = zeros(n,duration);   % Sucrose mass at the center of the cell
% ps1         = zeros(n,duration);   % Dynamic fluid pressure at the center of the cell
% vs1         = zeros(n,duration);   % Radial velocity at the center of the cell
% us1         = zeros(n-1,duration); % Axial velocity at the face of the cell
% nus1        = zeros(n,duration);   % Dynamic viscosity

%%% For saving
cs         = zeros(n,duration);   % Sucrose mass at the center of the cell
ps         = zeros(n,duration);   % Dynamic fluid pressure at the center of the cell
vs         = zeros(n,duration);   % Radial velocity at the center of the cell
us         = zeros(n-1,duration); % Axial velocity at the face of the cell
nus        = zeros(n,duration);   % Dynamic viscosity
cc         = zeros(n,duration);   % Sucrose concentration
pt         = zeros(n,duration);   % Total pressure (pd + ph)
cm         = zeros(1,duration);   % Total mass in the phloem 
dp         = zeros(1,duration);   % Dynamic pressure gradient along the phloem
psi_xylem  = zeros(1,duration);   % Xylem water potential
psi_phloem = zeros(1,duration);   % Phloem water potential
psi_ratio  = zeros(1,duration);   % Ratio
psi_diff   = zeros(1,duration);   % Difference
maxU       = zeros(1,duration);  % Maximum longitudinal velocity
cfl        = zeros(1,duration);  % Courant Number within each hour


% cmazen     = zeros(n,duration*nts);
% pmazen     = zeros(n,duration*nts);
% vmazen     = zeros(n,duration*nts);
% umazen     = zeros(n-1,duration*nts);
% numazen     = zeros(n,duration*nts);
% dpmazen     = zeros(1,duration*nts);

sinks = zeros(1,duration);
sinkr = zeros(1,duration);
source = zeros(1,duration);



fname = "./../results/paper/v5.xlsx"; % Filename to save

m     = 0;        % k = 0 constant viscosity, 
                  % k = 1 variable visvosity
k1    = 0;        % k1 = 0 constant sink, 
                  % k1 = 1 linear sink (highest at bottom), 
                  % k1 = 2 linear sink (highest at top),
                  % k1 = 3 2nd order polynomial sink (highest at bottom)
                  % k1 = 4 cts respiration and linear growth (highest at bottom)
                  % k1 = 5 cts respiration and linear growth (highest at top)
for j = 1:duration
    dx    = data.dbh(j);                 % diameter at breast height m
    L     = data.height(j);              % tree length m
    psi0  = - data.Psi_L(j);             % Leaf water potential Pa
    Sstem = ( data.Resp_s(j) ...
    + data.Grow_Resp_s(j) ...
    + data.Growth_s(j) )/(3600*L*pi*dx); % Stem sucrose sink g/(m s)
    Sstemr = ( data.Resp_s(j) ...
    + data.Grow_Resp_s(j) )/(3600*L*pi*dx); % Stem sucrose sink g/(m s)
    Sstemg = data.Growth_s(j)/(3600*L*pi*dx); % Stem sucrose sink g/(m s)
    Sroot = ( data.Resp_r(j) ...
    + data.Grow_Resp_r(j) ...
    + data.Growth_r(j) )/(3600*a*pi*dx); % Root sucrose sink g/(m s)
    

    
    Sleaf = ( Sstem*L + Sroot*a )/a;     % Sucrose source g/(m s)

    sinks(j)  = Sstem*L*3600*pi*dx ; 
    sinkr(j)  = Sroot*a*3600*pi*dx; 
    source(j) = sinks(j) + sinks(j);
    

    %%% Non-dimensional and scaling quantities
    c0 = psi0*Mw*pi*dx/(Rg*T);           % Sucrose scaling g/m2 
                                         % (assumed by balancing osmotic potential to leaf water potential)
                                         % it wont matter how you balance it
    cw = c0/(Mw*pi*dx);                  % mol/m3
    nu0 = 0.0017;
    % nu0 = exp(bb(1) + bb(2)*cw ...
    %     + bb(3)*(cw^2) + bb(4)*(cw^3) ...
    %     + bb(5)*(cw^4) );                % Viscosity of sucrose Pa*s
    os = Rg*T*c0/(Mw*pi*dx);             % Osmosis Pa
    es = a/L;                            % aspect ratio
    v0 = k*os;                           % Radial velocity scaling m/s
    u0 = v0/es;                          % Axial velocity scaling m/s
    p0 = L*nu0*u0/(a^2);                 % Pressure scaling Pa
    t0 = (a^2)/D;                        % diffusive timescale s
    G  = rho*g*L/os;                     % Gravity / osmosis
    X0 = psi0/os;                        % Xylem water potential / osmosis
    Mu = k*nu0*(L^2)/(a^3);              % Munch nb, axial / radial resistance
    Pe = a*v0/D;                         % Peclet number, advection / diffusion (radial)
    Ss = Sstem/(u0*c0);                  % Stem sink / advection
    Sl = Sleaf/(u0*c0);                  % Leaf source / advection
    Sr = Sroot/(u0*c0);                  % Root sink / advection

    Ss1 = Sstemr/(u0*c0);
    Ss2 = Sstemg/(u0*c0);

    %% Initial condition
    Psi = linspace(data.Psi_L(j)-rho*g*L,data.Psi_W(j),n)'./psi0; % Xylem water potential at the center of the cell
    if j==1
        [c(:,1),u(:,1),nu(:,1),v(:,1),p(:,1)] ...
            = initial(n,dz,Mu,G,X0,Pe,Sl,Ss,Sr,es,Psi, ...
            cw,bb,m,k1,Ss1,Ss2,dt,nts);
        % cm = dz*( sum(c(2:n-1,1)) + (1/2)*( c(1,1) + c(n,1) ) )

        % [cs1(:,j),us1(:,j),nus1(:,j),vs1(:,j),ps1(:,j)] ...
        %      = initial(n,dz,Mu,G,X0,Pe,Sl,Ss,Sr,es,Psi, ...
        %        cw,bb,m,k1,Ss1,Ss2,dt,nts);
        % % cm = dz*( sum(cs(2:n-1,1)) + (1/2)*( cs(1,1) + cs(n,1) ) )
           
    else
        c(:,1)  = co;
        u(:,1)  = uo;
        v(:,1)  = vo;
        p(:,1)  = po;
        nu(:,1) = nuo;
         
        % % Solve steady state
        % [cs1(:,j),us1(:,j),vs1(:,j),ps1(:,j),nus1(:,j),~] ...
        %     = iterate(n,dz,Mu,G,X0,Pe,Sl,Ss,Sr,es,...
        %       Psi,cs1(:,j-1),us1(:,j-1),vs1(:,j-1),ps1(:,j-1),nus1(:,j-1),...
        %       cw,bb,m,k1,0,dt,cs1(:,j-1),Ss1,Ss2,nts);
        % % cm = dz*( sum(cs(2:n-1,1)) + (1/2)*( cs(1,1) + cs(n,1) ) )
    end
    
    %% Solve in time
    for i = 2:nts
        [c(:,i),u(:,i),v(:,i),p(:,i),nu(:,i)] = ...
            iterate(n,dz,Mu,G,X0,Pe,Sl,Ss,Sr,es,...
            Psi,c(:,i-1),u(:,i-1),v(:,i-1),p(:,i-1),nu(:,i-1),...
            cw,bb,m,k1,1,dt,c(:,i-1),Ss1,Ss2,nts);

        % [c(:,i),u(:,i),v(:,i),p(:,i),nu(:,i)] ...
        %         = solve(n,dz,Mu,G,X0,Pe,Sl,Ss,Sr,es,Psi,...
        %         c(:,i-1),u(:,i-1),v(:,i-1),p(:,i-1),nu(:,i-1),...
        %         cw,bb,m,k1,dt,Ss1,Ss2,nts);

    %     cm = dz*( sum(c(2:n-1,i)) + (1/2)*( c(1,i) + c(n,i) ) )
    end
    
        
    % cs(:,((j-1)*nts+1):nts*j)         = c;
    % cmazen(:,((j-1)*nts+1):nts*j) = c;
    % pmazen(:,((j-1)*nts+1):nts*j) = p;
    % umazen(:,((j-1)*nts+1):nts*j) = u;
    % vmazen(:,((j-1)*nts+1):nts*j) = v;
    % numazen(:,((j-1)*nts+1):nts*j) = nu;
    % dpmazen(((j-1)*nts+1):nts*j) = p(1,:) - p(end,:);


    co              = c(:,end);
    uo              = u(:,end);
    vo              = v(:,end);
    po              = p(:,end);
    nuo             = nu(:,end);


    % cs(:,j)         = (1e-3).*cs1(:,j).*c0; 
    % ps(:,j)         = (1e-6).*ps1(:,j).*p0;
    % vs(:,j)         = vs1(:,j).*v0;
    % us(:,j)         = us1(:,j).*u0;
    % nus(:,j)        = nus1(:,end).*nu0;

    % make the variables dimensional
    cs(:,j)         = (1e-3).*co.*c0;   % Kg       
    ps(:,j)         = (1e-6).*po.*p0;   % MPa
    vs(:,j)         = vo.*v0;           % m/s
    us(:,j)         = uo.*u0;           % m/s
    nus(:,j)        = (1./nuo).*nu0;    % Pa s

    % variables for saving
    cc(:,j)         = (1e3).*cs(:,j)./(Mw*pi*dx);    % mol/m3
    cm(j)           = (1e-3)*pi*dx*L*dz*( sum(co(2:n-1))*c0 + (1/2)*( co(1)*c0 + co(n)*c0 ) ); % total mass in Kg
    % cm(j)           = (1e-3)*pi*dx*L*dz*( sum(cs1(2:n-1,j))*c0 + (1/2)*( cs1(1,j)*c0 + cs1(n,j)*c0 ) ); % total mass in Kg
    pt(:,j)         = (1e-6).*(p0.*po(:) + dz*L*rho*g.*((1:n)' - 1/2)); % total presusre in MPa (dynamic + hydrostatic)
    pp              = (1e-6).*(p0.*po(:) + dz*L*rho*g.*((1:n)' - 1/2) - os.*co(:)); % phloem water potential MPa
    % pt(:,j)         = (1e-6).*(p0.*ps1(:,j) + dz*L*rho*g.*((1:n)' - 1/2)); % total presusre in MPa (dynamic + hydrostatic)
    % pp              = (1e-6).*(p0.*ps1(:,j) + dz*L*rho*g.*((1:n)' - 1/2) - os.*cs1(:,j)); % phloem water potential MPa

    px              = (1e-6)*Psi.*psi0; % xylem water potential in MPa

    dp(j)           = (1e-6)*p0*(po(1) - po(end)); % pressure gradient in MPa
    % dp(j)           = (1e-6)*p0*(ps1(1,j) - ps1(end,j)); % pressure gradient in MPa

    psi_phloem(1,j) = dz*( sum(pp(2:n-1)) + (1/2)*( pp(1) + pp(n) ) ); % average phloem water potential MPa
    psi_xylem(1,j)  = dz*( sum(px(2:n-1)) + (1/2)*( px(1) + px(n) ) ); % average xylem water potential MPa
    psi_diff(1,j)   = psi_xylem(j) - psi_phloem(j); % average water potential disequilibrium Pa
    psi_ratio(1,j)  = psi_xylem(j)/psi_phloem(j);
    maxU(j)          = max(uo(:))*u0; 
    cfl(j)          = max(max(u))*dt/dz;
    % 
    % maxU(j)          = max(us1(:,j))*u0;
    % cfl(j)          = max(max(us1(:,1)))*dt/dz;


end

% PlotFigs(dt,duration,t0/3600,n,dz,L,cs,c0,vs,v0,us,u0,ps,p0,rho,g, ...
%     psi_phloem,psi_xylem,dp,psi_diff,psi_ratio)

% saveData(cs,us,vs,ps,nus,dp,psi_phloem,psi_xylem,pt,maxU,cc,cm,psi_diff,cfl,fname)
% exit

uavg1 = (1e6).*dp.*(a^2)./(3.*nus(1,:).*L);

% psipv = readmatrix("./../results/paper/" + "v5" + ".xlsx",'Sheet','PWP');
% cflv  = readmatrix("./../results/paper/" + "v5" + ".xlsx",'Sheet','maxU');
% ctv   = readmatrix("./../results/paper/" + "v5" + ".xlsx",'Sheet','TotalMass');
% uc    = readmatrix("./../results/paper/" + "v5" + ".xlsx",'Sheet','AxialV');
% vc    = readmatrix("./../results/paper/" + "v5" + ".xlsx",'Sheet','RadialV');
% pc    = readmatrix("./../results/paper/" + "v5" + ".xlsx",'Sheet','DynamicP');
% 
% 100.*(ps(:,end) - interp1(z1',pc(:,end),z2','spline','extrap'))./ps(:,end)
% 
% 
% z1 = (1/100).*(1/2:(100-1/2))
% z2 = dz.*((1/2):(n-1/2))
