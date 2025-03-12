function [fig1,fig2,fig3,fig4,fig5,fig6,fig7,fig8,fig9,fig10] =   ...
    PlotFigs(dt,period,t0,n,dz,L,cs,c0,vs,v0,us,u0,ps,p0,rho,g, ...
    psi_phloem,psi_xylem,dp,psi_diff,psi_ratio)

%%%% C at exit in time
fig1 = figure;
figure(fig1)
plot((1:period).*dt.*t0,cs(end,:).*c0,'Linewidth',10)
hold on
% plot((1:duration*nts).*dt.*t0,cs_v(end,:).*c0,'Linewidth',10)
% hold off
title('Concentration at the exit')
ax = findobj(fig1,'type','axes');
set(ax,'fontweight','bold','FontSize',55)
set([ax.XLabel],'string','{\boldmath$t (s)$}','Interpreter','latex')
set([ax.YLabel],'String','{\boldmath$C (g/m^2)$}','Interpreter','latex')

%%%% C at entrance in time
fig2 = figure;
figure(fig2)
plot((1:period).*dt.*t0,cs(1,:).*c0,'Linewidth',10)
hold on
% plot((1:duration*nts).*dt.*t0,cs_v(1,:).*c0,'Linewidth',10)
% hold off
title('Concentration at the entrance')
ax = findobj(fig2,'type','axes');
set(ax,'fontweight','bold','FontSize',55)
set([ax.XLabel],'string','{\boldmath$t (s)$}','Interpreter','latex')
set([ax.YLabel],'String','{\boldmath$C (g/m^2)$}','Interpreter','latex')

%%%% C profile at final time
fig3 = figure;
figure(fig3)
plot(((1:n)-1/2).*dz.*L,cs(:,end).*c0,'Linewidth',10)
hold on
% plot(((1:n)-1/2).*dz.*L,cs_v(:,end).*c0,'Linewidth',10)
% hold off
title('C profile at final time')
ax = findobj(fig3,'type','axes');
set(ax,'fontweight','bold','FontSize',55)
set([ax.XLabel],'string','{\boldmath$z (m)$}','Interpreter','latex')
set([ax.YLabel],'String','{\boldmath$C (g/m^2)$}','Interpreter','latex')

%%%% V profile at final time
fig4 = figure;
figure(fig4)
plot(((1:n) -1/2).*dz.*L,vs(:,end).*v0,'Linewidth',10)
title('V profile at final time')
ax = findobj(fig4,'type','axes');
set(ax,'fontweight','bold','FontSize',55)
set([ax.XLabel],'string','{\boldmath$z (m)$}','Interpreter','latex')
set([ax.YLabel],'String','{\boldmath$v (m/s)$}','Interpreter','latex')
 
%%%% U profile at final time
fig5 = figure;
figure(fig5)
plot((1:n-1).*dz.*L,us(:,end).*u0,'Linewidth',10)
title('Axial velocity profile at t = 3600 s')
ax = findobj(fig5,'type','axes');
set(ax,'fontweight','bold','FontSize',55)
set([ax.XLabel],'string','{\boldmath$z (m)$}','Interpreter','latex')
set([ax.YLabel],'String','{\boldmath$u (m/s)$}','Interpreter','latex')


%%%% Phloem and Xylem water potential (and pressure)
fig6 = figure;
figure(fig6)
plot((1:period),psi_phloem,'Linewidth',10)
hold on
plot((1:period),psi_xylem,'Linewidth',10)
hold off
title('Water Potential')
ax = findobj(fig6,'type','axes');
set(ax,'fontweight','bold','FontSize',55)
set([ax.XLabel],'string','{\boldmath$t (m)$}','Interpreter','latex')
set([ax.YLabel],'String','{\boldmath$\Psi$}','Interpreter','latex')

%%%% Pressure gradient in time
fig7 = figure;
figure(fig7)
plot((1:period),dp,'Linewidth',10)
title('Water Potential')
ax = findobj(fig7,'type','axes');
set(ax,'fontweight','bold','FontSize',55)
set([ax.XLabel],'string','{\boldmath$t (hr)$}','Interpreter','latex')
set([ax.YLabel],'String','{\boldmath$\Delta p / L$}','Interpreter','latex')

%%%% Water potential difference in time
fig8 = figure;
figure(fig8)
plot((1:period),psi_diff,'Linewidth',10)
hold on
title('Water Potential')
ax = findobj(fig8,'type','axes');
set(ax,'fontweight','bold','FontSize',55)
set([ax.XLabel],'string','{\boldmath$t (hr)$}','Interpreter','latex')
set([ax.YLabel],'String','{\boldmath$\Psi_x - \Psi_p$}','Interpreter','latex')

%%%% Water potential ratio in time
fig9 = figure;
figure(fig9)
plot((1:period),psi_ratio,'Linewidth',10)
hold on
title('Water Potential')
ax = findobj(fig9,'type','axes');
set(ax,'fontweight','bold','FontSize',55)
set([ax.XLabel],'string','{\boldmath$t (hr)$}','Interpreter','latex')
set([ax.YLabel],'String','{\boldmath$Ratio$}','Interpreter','latex')

%%%% P profile at final time
fig10 = figure;
figure(fig10)
plot(((1:n) -1/2).*dz.*L,ps(:,end).*p0 + dz*L*rho*g.*((1:n)' -1/2),'Linewidth',10)
title('P profile at final time')
ax = findobj(fig10,'type','axes');
set(ax,'fontweight','bold','FontSize',55)
set([ax.XLabel],'string','{\boldmath$z (m)$}','Interpreter','latex')
set([ax.YLabel],'String','{\boldmath$p (Pa)$}','Interpreter','latex')

end