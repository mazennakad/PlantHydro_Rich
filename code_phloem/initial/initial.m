function [c,u,nu,v,p] = initial(n,dz,Mu,G,X0,Pe,Sl,Ss,Sr,es,Psi,...
            cw,bb,m,k1,Ss1,Ss2,dt,nts)

[co,uo,vo,po,nuo] = guess(n,dz,G,Mu,X0,Psi,cw,bb,m);
[c,u,v,p,nu,~]    = iterate(n,dz,Mu,G,X0,Pe,Sl,Ss,Sr,es,...
                            Psi,co,uo,vo,po,nuo,cw,bb,m,k1,0,dt,zeros(n,1),Ss1,Ss2,nts);

end