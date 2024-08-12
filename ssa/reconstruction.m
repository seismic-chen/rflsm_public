function [d] = reconstruction(din,T,dt,fmin,fmax,rank_p,alpha,n_iter)

aux = din;

for k = 1:n_iter

    if k == n_iter
 	   alpha=0;
    end

    out = fx_ssa(aux,dt,rank_p,fmin,fmax);

    d = alpha*din + (1-alpha*T).*out;

    aux = d;

end

