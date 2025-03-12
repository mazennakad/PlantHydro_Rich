function [ck,uk,vk,pk,nuk,iteration] = iterate(n,dz,Mu,G,X0,Pe,Sl,Ss,Sr,es,...
        Psi,co,uo,vo,po,nuo,cw,bb,k,k1,s,dt,cp,Ss1,Ss2,nts)
check_C = true; % for stopping the iteration
iteration = 1;  % for the number of iterations
while check_C
    % Set the new variables
    if iteration ~= 1
        co = ck;
        uo = uk;
        vo  = vk;
        po  = pk;
        nuo = nuk;
    end

    [vector_S,Sk] = newton(n,dz,Mu,G,X0,Pe,Sl,Ss,Sr,es,...
        Psi,co,uo,vo,po,nuo,cw,bb,k,k1,s,dt,cp,Ss1,Ss2,nts);
    check_C = check_convergence(vector_S,Sk,length(Sk));
    pk  = Sk(1:n);
    uk  = Sk(n+1:2*n-1);
    vk  = Sk(2*n:3*n-1);
    ck  = Sk(3*n:4*n-1);
    nuk = Sk(4*n:end);
    if iteration > 50
        break;
    end
    iteration = iteration + 1; % to set a maximum number of iterations
end

end
