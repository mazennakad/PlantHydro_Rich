function saveData(c,u,v,p,nu,dp,psip,psix,pt,maxU,cc,cm,psi_diff,cfl,filename)
% c sugar mass g/m2
% u axial velocity m/s
% v radial velocity m/s
% p dynamic pressure Pa
% nu dynamic viscosity Pa/(m s)
% dp dyanmic pressure gradient Pa/m
% psip and psix average phloem and xylem water potential Pa
% pt total pressure Pa
% maxU is the maximum velocity at each hour
% cc sugar concentration mol/m3
% cm total mass g
% psi_diff is the difference in water potential
% cfl is the courant number

writematrix(c,filename,'Sheet','Mass')
writematrix(u,filename,'Sheet','AxialV')
writematrix(v,filename,'Sheet','RadialV')
writematrix(p,filename,'Sheet','DynamicP')
writematrix(nu,filename,'Sheet','DynamicV')
writematrix(dp,filename,'Sheet','PressureG')
writematrix(psip,filename,'Sheet','PWP')
writematrix(psix,filename,'Sheet','XWP')
writematrix(pt,filename,'Sheet','TotalP')
writematrix(maxU,filename,'Sheet','MaxU')
writematrix(cc,filename,'Sheet','Concentration')
writematrix(cm,filename,'Sheet','TotalMass')
writematrix(psi_diff,filename,'Sheet','DWP')
writematrix(cfl,filename,'Sheet','cfl')

end