more off;
addpath("/auto/jlischner/Dropbox/EG2d/code/fullfreq");
angstrom = 1/0.5291;
sec = 1/2.418884326505e-17; %# wikipedia: atomic units
cm = 10^8*angstrom;

load ffoutput;

#global nv = 1;
#global m  = 0.067;
#global kF = sqrt( n_cm/nv/cm^2 *2*pi );
#global EF = kF^2/2/m;

nus = ws/EF;
wss = [min(nus)+0.01 : 0.01 : max(nus)-0.01]'*EF;
Nfreq = length(wss);
Nk = length(kys);
Ass = zeros(Nfreq,Nk);
Sigss = zeros(Nfreq,Nk);

%# get dmu:
[a,b] = min(abs(kys-1));
if( abs(a) < 0.000001);
  ikf = b;
else;
  printf("could not find kF in kys!");
  exit;
endif;

[a,b] = min(abs(nus));
if( abs(a) < 0.001 );
  dmu = Sigs(b,ikf);
  printf("dmu = %f meV \n",dmu*27.21*1000);
else;
  printf("w=0 not in frequency grid! \n");
  exit;
endif;

%# problem in GW+C at k close to kF
Nkmax = 1;

for ik = 1:Nkmax;
  
  printf("doing ik=%d of %d states \n",ik,Nkmax);
  k = kys(ik)*kF;
  xi = k^2/2/m - EF;
  Sig = Sigs(:,ik);
  ReSig = interp1(ws,real(Sig),wss);
  ImSig = interp1(ws,imag(Sig),wss);

  [a,b] = max(abs(ImSig));
  shift = ImSig(b)/100;
  ImSig += shift;

  Sig = ReSig + I*ImSig;
  A   = 1/pi * abs(imag(Sig)) ./ ( (wss-xi-real(Sig-dmu)).^2 + imag(Sig).^2 );
  
  Ass(:,ik) = A;
  kAss(:,ik) = k*A;
  Sigss(:,ik) = Sig;
  
  dat = 27.21*[wss,ReSig,ImSig,ImSig];
  save gwc_input.dat dat;
  [wsc,Ac] = get_gwc(xi*27.21,dmu*27.21);
  
  %# normalize GW+C spectral function to same value as GW
  Ac = real(Ac);
  normAc = trapz(wsc,Ac);
  normA  = trapz(wss,A);
  Ac *= normA/normAc;

  Ac *= 27.21; %# convert to a.u.
  Acs(:,ik)  = Ac;
  kAcs(:,ik) = k*Ac;

end;

dk = kF/Nk;
dos  = sum(kAss,2) * dk/pi;
dosc = sum(kAcs,2) * dk/pi;

wss *= 27.21*1000; %# convert to meV
wsc *= 1000;

save gwc_output wss wsc dos dosc Acs Ass;

more on;