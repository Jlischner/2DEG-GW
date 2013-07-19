function [out,chi0,chi] = diel_eta(q,w);
    
    global kF;
    global vF;
    global EF;
    global eps_back;
    global m;
    global Vc;
    #Vq = 2*pi./q/eps_back;
    N0  = m/pi;
    #tol = 0.8;
    #eta = 2*EF;
    
    %# from Vignale book
    qb = q/kF;
    w = abs(real(w)) + I*abs(imag(w));
    vp = w./(q*vF) + q./(2*kF);
    vm = w./(q*vF) - q./(2*kF);
    rvp = real(vp);
    rvm = real(vm);

    Rchi = sign(rvm) .* (rvm.^2 > 1) .*sqrt(vm.^2-1) - sign(rvp) .* ...
           (rvp.^2 >1 ) .* sqrt(vp.^2 - 1);
    Rchi = - N0 * (1 + Rchi ./qb);
    
    Ichi = (1>rvm.^2) .* sqrt(1-vm.^2) - (1>rvp.^2) .* sqrt(1-vp.^2);
    Ichi = - N0*Ichi./qb;
    
    chi0 = Rchi + I*Ichi;
    
    out = 1 - Vc .* chi0;
#    out(real(w)<0) = conj(out(real(w)<0));
#    out( abs(real(out)) < tol) += I*eta;
    chi = (1./out -1)./Vc;
endfunction;