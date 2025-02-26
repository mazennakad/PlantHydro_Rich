clear all
close all

path   = "./../results/paper/";
% datac   = ["c0","c1","c2","c4","c5"];
datac   = ["c0","c1","c2","c4","c5","rs"];
datav   = ["v0","v1","v2","v4","v5"];
color  = [[0 0 1]; [0.5 0.0 0.5]; [1 0 1]; ...
         [0.5 0.5 0.5]; [0.5 0 0]; [0 0.5 0.5]];
n      = 100;
dz     = 1/n;
a     = 1e-5;                        % phloem thickness m
Mw    = 342.2965;                    % Molecular mass of sucrose g/mol



file   = readtable("/Users/mazennakad/Desktop/phloem/phloem_biomeE" + ...
        "/output/BCI_hydro_Cohort_hourly.csv");
L      = file.height(1);
dx    = file.dbh(1:60);                 % diameter at breast height m


t      = 1:24;   
figC   = figure;
figU   = figure;
% figV   = figure;
% figP   = figure;
figPg  = figure;
% figPsi = figure;
% figCfl = figure;
% figCt  = figure;

for i = 1:6
    cc    = readmatrix(path + datac(i) + ".xlsx",'Sheet','Mass');
    uc    = readmatrix(path + datac(i) + ".xlsx",'Sheet','AxialV');
    vc    = readmatrix(path + datac(i) + ".xlsx",'Sheet','RadialV');
    pc    = readmatrix(path + datac(i) + ".xlsx",'Sheet','DynamicP');
    nuc   = readmatrix(path + datac(i) + ".xlsx",'Sheet','DynamicV');
    dpc   = readmatrix(path + datac(i) + ".xlsx",'Sheet','PressureG');
    psipc = readmatrix(path + datac(i) + ".xlsx",'Sheet','PWP');
    psix  = readmatrix(path + datac(i) + ".xlsx",'Sheet','XWP');
    ptc   = readmatrix(path + datac(i) + ".xlsx",'Sheet','TotalP');
    cflc  = readmatrix(path + datac(i) + ".xlsx",'Sheet','maxU');
    ctc   = readmatrix(path + datac(i) + ".xlsx",'Sheet','TotalMass');
    
    % cv    = readmatrix(path + datav(i) + ".xlsx",'Sheet','Mass');
    % uv    = readmatrix(path + datav(i) + ".xlsx",'Sheet','AxialV');
    % vv    = readmatrix(path + datav(i) + ".xlsx",'Sheet','RadialV');
    % pv    = readmatrix(path + datav(i) + ".xlsx",'Sheet','DynamicP');
    % nuv   = readmatrix(path + datav(i) + ".xlsx",'Sheet','DynamicV');
    % dpv   = readmatrix(path + datav(i) + ".xlsx",'Sheet','PressureG');
    % psipv = readmatrix(path + datav(i) + ".xlsx",'Sheet','PWP');
    % ptv   = readmatrix(path + datav(i) + ".xlsx",'Sheet','TotalP');
    % cflv  = readmatrix(path + datav(i) + ".xlsx",'Sheet','maxU');
    % ctv   = readmatrix(path + datav(i) + ".xlsx",'Sheet','TotalMass');

    % cv    = readmatrix(path1 + datac(i) + ".xlsx",'Sheet','Mass');
    % uv    = readmatrix(path1 + datac(i) + ".xlsx",'Sheet','AxialV');
    % vv    = readmatrix(path1 + datac(i) + ".xlsx",'Sheet','RadialV');
    % pv    = readmatrix(path1 + datac(i) + ".xlsx",'Sheet','DynamicP');
    % nuv   = readmatrix(path1 + datac(i) + ".xlsx",'Sheet','DynamicV');
    % dpv   = readmatrix(path1 + datac(i) + ".xlsx",'Sheet','PressureG');
    % psipv = readmatrix(path1 + datac(i) + ".xlsx",'Sheet','PWP');
    % ptv   = readmatrix(path1 + datac(i) + ".xlsx",'Sheet','TotalP');
    % cflv  = readmatrix(path1 + datac(i) + ".xlsx",'Sheet','maxU');
    % ctv   = readmatrix(path1 + datac(i) + ".xlsx",'Sheet','TotalMass');
    % 
    % cv1    = readmatrix(path2 + datac(i) + ".xlsx",'Sheet','Mass');
    % uv1    = readmatrix(path2 + datac(i) + ".xlsx",'Sheet','AxialV');
    % vv1    = readmatrix(path2 + datac(i) + ".xlsx",'Sheet','RadialV');
    % pv1    = readmatrix(path2 + datac(i) + ".xlsx",'Sheet','DynamicP');
    % nuv1   = readmatrix(path2 + datac(i) + ".xlsx",'Sheet','DynamicV');
    % dpv1   = readmatrix(path2 + datac(i) + ".xlsx",'Sheet','PressureG');
    % psipv1 = readmatrix(path2 + datac(i) + ".xlsx",'Sheet','PWP');
    % ptv1   = readmatrix(path2 + datac(i) + ".xlsx",'Sheet','TotalP');
    % cflv1  = readmatrix(path2 + datac(i) + ".xlsx",'Sheet','maxU');
    % ctv1   = readmatrix(path2 + datac(i) + ".xlsx",'Sheet','TotalMass');
    
    cavg = dz.*(sum(cc(2:n-1,:)) + (1/2).*( cc(1,:) + cc(n,:) ))./(pi.*dx(1:24)'.*Mw*1e-3);
    uavg = dz.*(sum(uc(2:n-2,:)) + (1/2).*( uc(1,:) + uc(n-1,:) ));
    % kavg = uavg./dpc;
    % uavg1 = (1e6).*dpc.*(a^2)./(3.*nuc(1v,:).*L);
    figure(figC)
    plot(t,cavg,'-','Color',color(i,:),'Linewidth',3)
    hold on
    % plot(((1:n)-1/2).*dz.*L,cc(:,24),'-','Color',color(i,:),'Linewidth',3)
    % hold on
    % plot(((1:n)-1/2).*dz.*L,cv(:,24),'--','Color',color(i,:),'Linewidth',3)%,'HandleVisibility','off')
    % hold on
    % plot(((1:n)-1/2).*dz.*L,cc(:,18),'-','Color',color(1,:),'Linewidth',3)
    % hold on
    % plot(((1:n)-1/2).*dz.*L,cc(:,12),'-','Color',color(2,:),'Linewidth',3)
    % hold on
    % plot(((1:n)-1/2).*dz.*L,cc(:,6),'-','Color',color(3,:),'Linewidth',3)
    % hold on
    

    figure(figU)
    plot(t,uavg,'-','Color',color(i,:),'Linewidth',3)
    hold on
    % % plot(((1:n-1)).*dz.*L,uc(:,24),'-','Color',color(i,:),'Linewidth',3)
    % % hold on
    % plot(((1:n-1)).*dz.*L,uv(:,24),'--','Color',color(i,:),'Linewidth',3)%,'HandleVisibility','off')
    % hold on
    % plot(((1:n-1)).*dz.*L,uc(:,18),'-','Color',color(1,:),'Linewidth',3)
    % hold on
    % plot(((1:n-1)).*dz.*L,uc(:,12),'-','Color',color(2,:),'Linewidth',3)
    % hold on
    % plot(((1:n-1)).*dz.*L,uc(:,6),'-','Color',color(3,:),'Linewidth',3)
    % hold on


    % 
    % 
    % figure(figV)
    % plot(((1:n)-1/2).*dz.*L,vc(:,24),'-','Color',color(i,:),'Linewidth',3)
    % hold on
    % plot(((1:n)-1/2).*dz.*L,vv(:,24),'--','Color',color(i,:),'Linewidth',3)%,'HandleVisibility','off')
    % hold on
    % plot(((1:n)-1/2).*dz.*L,vc(:,18),'-','Color',color(1,:),'Linewidth',3)
    % hold on
    % plot(((1:n)-1/2).*dz.*L,vc(:,12),'-','Color',color(2,:),'Linewidth',3)
    % hold on
    % plot(((1:n)-1/2).*dz.*L,vc(:,6),'-','Color',color(3,:),'Linewidth',3)
    % hold on

 

    % figure(figP)
    % plot(((1:n)-1/2).*dz.*L,ptc(:,24),'-','Color',color(i,:),'Linewidth',3)
    % hold on
    % plot(((1:n)-1/2).*dz.*L,ptv(:,24),'--','Color',color(i,:),'Linewidth',3)%,'HandleVisibility','off')
    % hold on
    % plot(((1:n)-1/2).*dz.*L,ptc(:,18),'-','Color',color(1,:),'Linewidth',3)
    % hold on
    % plot(((1:n)-1/2).*dz.*L,ptc(:,12),'-','Color',color(2,:),'Linewidth',3)
    % hold on
    % plot(((1:n)-1/2).*dz.*L,ptc(:,6),'-','Color',color(3,:),'Linewidth',3)
    % hold on

 

    % figure(figPsi)
    % % % plot(t,psix(t),'-','Color',[0 0 0],'Linewidth',3,'HandleVisibility','off')
    % % % hold on
    % % plot(t,psipc(t),'-','Color',color(i,:),'Linewidth',3)
    % % hold on
    % plot(t,psipv(t),'--','Color',color(i,:),'Linewidth',3)%,'HandleVisibility','off')
    % hold on
    

    % figure(figCt)
    % % plot(t,ctc(t),'-','Color',color(i,:),'Linewidth',3)
    % % hold on
    % plot(t,ctv(t),'--','Color',color(i,:),'Linewidth',3)%,'HandleVisibility','off')
    % hold on


    figure(figPg)
    plot(t,dpc(t),'-','Color',color(i,:),'Linewidth',3)
    hold on
    % plot(t,dpv(t),'--','Color',color(i,:),'Linewidth',3)%,'HandleVisibility','off')
    % hold on
    % 


    % figure(figCfl)
    % % plot(t,cflc(t),'-','Color',color(i,:),'Linewidth',3)
    % % hold on
    % plot(t,cflv(t),'--','Color',color(i,:),'Linewidth',3)%,'HandleVisibility','off')
    % hold on

end

figure(figC)
ax = findobj(figC,'type','axes');
set(ax,'fontweight','bold','FontSize',55)
set([ax.XLabel],'string','{\boldmath$t (hr)$}','Interpreter','latex')
set([ax.YLabel],'String','{\boldmath$\bar{c} (mol/m^3)$}','Interpreter','latex')
legend('CS','LSb','LSt','LGSb','LGSt','RS','location','northeast')
figure(figU)
ax = findobj(figU,'type','axes');
set(ax,'fontweight','bold','FontSize',55)
set([ax.XLabel],'string','{\boldmath$t (hr)$}','Interpreter','latex')
set([ax.YLabel],'String','{\boldmath$\bar{u} (m/s)$}','Interpreter','latex')
legend('CS','LSb','LSt','LGSb','LGSt','RS','location','northeast')
figure(figPg)
ax = findobj(figPg,'type','axes');
set(ax,'fontweight','bold','FontSize',55)
set([ax.XLabel],'string','{\boldmath$t (hr)$}','Interpreter','latex')
set([ax.YLabel],'String','{\boldmath$\Delta p (MPa)$}','Interpreter','latex')
legend('CS','LSb','LSt','LGSb','LGSt','RS','location','northeast')



% 
% 
% figure(figC)
% title('C profile')
% ax = findobj(figC,'type','axes');
% set(ax,'fontweight','bold','FontSize',55)
% set([ax.XLabel],'string','{\boldmath$z (m)$}','Interpreter','latex')
% set([ax.YLabel],'String','{\boldmath$C (Kg/m^2)$}','Interpreter','latex')
% legend('CS','LSb','LSt','LGSb','LGSt','location','northeast')
% % legend('SS','TRI','TRS','location','northeast')
% % legend('CS','LSt','LGSt','location','northeast')
% % legend('24','18','12','6','location','northeast')
% 
% 
% 
% 
% figure(figPsi)
% title('Water Potential')
% ax = findobj(figPsi,'type','axes');
% set(ax,'fontweight','bold','FontSize',55)
% set([ax.XLabel],'string','{\boldmath$t (hr)$}','Interpreter','latex')
% set([ax.YLabel],'String','{\boldmath$\Psi (MPa)$}','Interpreter','latex')
% legend('CS','LSb','LSt','LGSb','LGSt','location','northwest')
% % legend('SS','TRI','TRS','location','northeast')
% % legend('CS','LSt','LGSt','location','northeast')
% 
% 
% figure(figCt)
% title('Total Mass')
% ax = findobj(figCt,'type','axes');
% set(ax,'fontweight','bold','FontSize',55)
% set([ax.XLabel],'string','{\boldmath$t (hr)$}','Interpreter','latex')
% set([ax.YLabel],'String','{\boldmath$C_t (Kg)$}','Interpreter','latex')
% % legend('SS','TRI','TRS','location','southwest')
% legend('CS','LSb','LSt','LGSb','LGSt','location','northwest')
% % legend('CS','LSt','LGSt','location','northeast')
% 
% 
% figure(figV)
% yline(0,'--','Color',[0 0 0],'Linewidth',2,'HandleVisibility','off')
% title('Radial velocity profile')
% ax = findobj(figV,'type','axes');
% set(ax,'fontweight','bold','FontSize',55)
% set([ax.XLabel],'string','{\boldmath$z (m)$}','Interpreter','latex')
% set([ax.YLabel],'String','{\boldmath$v (m/s)$}','Interpreter','latex')
% legend('CS','LSb','LSt','LGSb','LGSt','location','southwest')
% % legend('SS','TRI','TRS','location','northeast')
% % legend('CS','LSt','LGSt','location','northeast')
% % legend('24','18','12','6','location','northeast')
% 
% 
% 
% figure(figU)
% title('Axial velocity profile')
% ax = findobj(figU,'type','axes');
% set(ax,'fontweight','bold','FontSize',55)
% set([ax.XLabel],'string','{\boldmath$z (m)$}','Interpreter','latex')
% set([ax.YLabel],'String','{\boldmath$u (m/s)$}','Interpreter','latex')
% legend('CS','LSb','LSt','LGSb','LGSt','location','southwest')
% % legend('SS','TRI','TRS','location','northeast')
% % legend('CS','LSt','LGSt','location','northeast')
% % legend('24','18','12','6','location','northeast')
% 
% 
% figure(figP)
% title('Pressure profile')
% ax = findobj(figP,'type','axes');
% set(ax,'fontweight','bold','FontSize',55)
% set([ax.XLabel],'string','{\boldmath$z (m)$}','Interpreter','latex')
% set([ax.YLabel],'String','{\boldmath$p (MPa)$}','Interpreter','latex')
% legend('CS','LSb','LSt','LGSb','LGSt','location','northwest')
% % legend('SS','TRI','TRS','location','northeast')
% % legend('CS','LSt','LGSt','location','northeast')
% % legend('24','18','12','6','location','northeast')
% 
% 
% figure(figPg)
% title('Pressure Gradient')
% ax = findobj(figPg,'type','axes');
% set(ax,'fontweight','bold','FontSize',55)
% set([ax.XLabel],'string','{\boldmath$t (hr)$}','Interpreter','latex')
% set([ax.YLabel],'String','{\boldmath$\Delta p (MPa)$}','Interpreter','latex')
% legend('CS','LSb','LSt','LGSb','LGSt','location','northwest')
% % legend('SS','TRI','TRS','location','northeast')
% % legend('CS','LSt','LGSt','location','northeast')
% 
% 
% figure(figCfl)
% title('Maximum velocity')
% ax = findobj(figCfl,'type','axes');
% set(ax,'fontweight','bold','FontSize',55)
% set([ax.XLabel],'string','{\boldmath$t (hr)$}','Interpreter','latex')
% set([ax.YLabel],'String','{\boldmath$U_{max} (m/s)$}','Interpreter','latex')
% % legend('SS','TRI','TRS','location','southwest')
% legend('CS','LSb','LSt','LGSb','LGSt','location','southwest')
% % legend('CS','LSt','LGSt','location','northeast')











% cc    = readmatrix(path + "c1.xlsx",'Sheet','Mass');
% uc    = readmatrix(path + "c1.xlsx",'Sheet','AxialV');
% vc    = readmatrix(path + "c1.xlsx",'Sheet','RadialV');
% pc    = readmatrix(path + "c1.xlsx",'Sheet','DynamicP');
% nuc   = readmatrix(path + "c1.xlsx",'Sheet','DynamicV');
% dpc   = readmatrix(path + "c1.xlsx",'Sheet','PressureG');
% psipc = readmatrix(path + "c1.xlsx",'Sheet','PWP');
% psix  = readmatrix(path + "c1.xlsx",'Sheet','XWP');
% ptc   = readmatrix(path + "c1.xlsx",'Sheet','TotalP');
% cflc  = readmatrix(path + "c1.xlsx",'Sheet','CourantNb');
% ctc   = readmatrix(path + "c1.xlsx",'Sheet','TotalMass');

% cv    = readmatrix(path + "v1.xlsx",'Sheet','Mass');
% uv    = readmatrix(path + "v1.xlsx",'Sheet','AxialV');
% vv    = readmatrix(path + "v1.xlsx",'Sheet','RadialV');
% pv    = readmatrix(path + "v1.xlsx",'Sheet','DynamicP');
% nuv   = readmatrix(path + "v1.xlsx",'Sheet','DynamicV');
% dpv   = readmatrix(path + "v1.xlsx",'Sheet','PressureG');
% psipv = readmatrix(path + "v1.xlsx",'Sheet','PWP');
% ptv   = readmatrix(path + "v1.xlsx",'Sheet','TotalP');
% cflv  = readmatrix(path + "v1.xlsx",'Sheet','CourantNb');
% ctv   = readmatrix(path + "v1.xlsx",'Sheet','TotalMass');

% figure(figC)
% a1 = axes('Position',[.65 .65 .25 .25])
% box on
% plot(((1:n)-1/2).*dz.*L,cc(:,end),'-','Color',color(i,:),'Linewidth',3)
% hold on
% plot(((1:n)-1/2).*dz.*L,cv(:,end),'--','Color',color(i,:),'Linewidth',3,'HandleVisibility','off')
% hold off
% legend('LSb','location','northwest')
% set(a1,'fontweight','bold','FontSize',55)
% 
% figure(figU)
% a1 = axes('Position',[.35 .65 .25 .25])
% box on
% plot(((1:n-1)).*dz.*L,uc(:,end),'-','Color',color(i,:),'Linewidth',3)
% hold on
% plot(((1:n-1)).*dz.*L,uv(:,end),'--','Color',color(i,:),'Linewidth',3,'HandleVisibility','off')
% hold off
% legend('LSb','location','northwest')
% set(a1,'fontweight','bold','FontSize',55)
% 
% 
% figure(figV)
% a1 = axes('Position',[.35 .65 .25 .25])
% box on
% plot(((1:n)-1/2).*dz.*L,vc(:,end),'-','Color',color(i,:),'Linewidth',3)
% hold on
% plot(((1:n)-1/2).*dz.*L,vv(:,end),'--','Color',color(i,:),'Linewidth',3,'HandleVisibility','off')
% hold off
% legend('LSb','location','northwest')
% set(a1,'fontweight','bold','FontSize',55)
% 
% figure(figP)
% a1 = axes('Position',[.65 .65 .25 .25])
% box on
% plot(((1:n)-1/2).*dz.*L,ptc(:,end),'-','Color',color(i,:),'Linewidth',3)
% hold on
% plot(((1:n)-1/2).*dz.*L,ptv(:,end),'--','Color',color(i,:),'Linewidth',3,'HandleVisibility','off')
% hold off
% legend('LSb','location','northwest')
% set(a1,'fontweight','bold','FontSize',55)
% 
% figure(figPsi)
% a1 = axes('Position',[.65 .65 .25 .25])
% box on
% plot(t,psix,'-','Color',[0 0 0],'Linewidth',3,'HandleVisibility','off')
% hold on
% plot(t,psipc,'-','Color',color(i,:),'Linewidth',3)
% hold on
% plot(t,psipv,'--','Color',color(i,:),'Linewidth',3,'HandleVisibility','off')
% hold off
% legend('LSb','location','northwest')
% set(a1,'fontweight','bold','FontSize',55)
% 
% figure(figCt)
% a1 = axes('Position',[.45 .65 .25 .25])
% box on
% plot(t,ctc,'-','Color',color(i,:),'Linewidth',3)
% hold on
% plot(t,ctv,'--','Color',color(i,:),'Linewidth',3,'HandleVisibility','off')
% hold off
% legend('LSb','location','northwest')
% set(a1,'fontweight','bold','FontSize',55)
% 
% figure(figPg)
% a1 = axes('Position',[.65 .65 .25 .25])
% box on
% plot(t,dpc,'-','Color',color(i,:),'Linewidth',3)
% hold on
% plot(t,dpv,'--','Color',color(i,:),'Linewidth',3,'HandleVisibility','off')
% hold off
% legend('LSb','location','northwest')
% set(a1,'fontweight','bold','FontSize',55)
% 
% figure(figCfl)
% a1 = axes('Position',[.65 .65 .25 .25])
% box on
% plot(t,cflc,'-','Color',color(i,:),'Linewidth',3)
% hold on
% plot(t,cflv,'--','Color',color(i,:),'Linewidth',3,'HandleVisibility','off')
% hold off
% legend('LSb','location','northwest')
% set(a1,'fontweight','bold','FontSize',55)
% 
% 
% 
