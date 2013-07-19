more off;
addpath("/auto/jlischner/Dropbox/EG2d/code/fullfreq");
angstrom = 1/0.5291;
sec = 1/2.418884326505e-17; %# wikipedia: atomic units
cm = 10^8*angstrom;

%# read input
input;
global nv = 1;
global kF = sqrt( n_cm/nv/cm^2 *2*pi );
global EF = kF^2/2/m;
global vF = kF/m;
global n = kF^2/2/pi*nv;
global eta = etaEF*EF;
rs = sqrt(2)/kF; %# from Polini/Vignale
kc = Lambda*kF; %# cutoff for correlation self energy

%# good convergence parameters:
%# dphi = 0.01 and du=0.001
phis = [0 : 0.01 : 2*pi]';
us   = [0.01 : 0.001 : 1]';
dphi = phis(2)-phis(1);
qs   = kc*us.^2; %# dq = 2*kc*us*dus
Nq   = length(qs);

%# finite well form factor (Das Sarma, Eq.31):
qd = qs*d;
fq = 3/8*qd + pi^2./qd - 4*pi^4*(1-exp(-qd))./(qd.^2 .*(qd.^2 + 4*pi^2) );
fq .*= 8./(qd.^2 + 4*pi^2);
if( abs(d) < 0.000001);
  fq = 1;
endif;

%# form factor due to metallic electrode (image charge model)
f_el = (1-exp(-2*d_el*qs));
if( abs(d_el) < 0.000001);
  f_el = 1;
endif;

%# set up screened 2D Coulomb interaction
global Vc = 2*pi./abs(qs) ./eps_back .*fq .* f_el;
qVc = qs .* Vc;

%# calculate self energies
ws = nus*EF;
Nfreq = length(ws);
Nk  = length(kys);
Sigs = zeros(Nfreq,Nk);
As   = zeros(Nfreq,Nk);

for ik = 1:Nk;

  k  = [0 kys(ik)*kF];
  SigXw = [];
  SigCw = [];

  for iw = 1:Nfreq;

    printf("state: %d of %d - freq: %d of %d \n",ik,Nk,iw,Nfreq);
    w = ws(iw);
    
    SigX = 0.;
    SigC = 0.;
    
    for phi = phis';
      
      kmq = -qs*[cos(phi) sin(phi)] + ones(Nq,1)*k;
      abs_kmq = sqrt(sum(kmq.^2,2));
      xi_kmq  = abs_kmq.^2/2/m - EF;
      SigX += trapz(us, us .* qVc .* (xi_kmq<0) );
            
      wts = sign(xi_kmq) .*(w-xi_kmq);
      [eps2d,chi0,chi] = diel_eta(qs, wts + I*eta);
      epsI = 1./eps2d -1;
      SigC += trapz(us, us .* qVc .* epsI .* (wts>0) );
      
    endfor;%# end phis loop
    
    SigX *= -kc*dphi*2/(2*pi)^2; 
    SigC *=  kc*dphi*2/(2*pi)^2;             
    
    SigXw = [SigXw; SigX];
    SigCw = [SigCw; SigC];
    
  endfor;%# end nus loop
  
  ImSig = imag(SigCw);
  ReSig = zeros(size(ImSig));
  Nfreq = length(ws);
  %# Kramer-Kronig: Vignale book (eq.6.80)
  for ii = 1:Nfreq;
    ReSig(ii) += trapz(ws, ImSig./(ws-ws(ii)+I*eta));
  endfor;
  ReSig = real(ReSig)/pi + SigXw;

  Sig = ReSig+I*ImSig;
  Sigs(:,ik) = Sig;

endfor;%# end k loop

ws = nus*EF;
save ffoutput ws Sigs kys m n_cm eps_back nv kF EF;
more on;
