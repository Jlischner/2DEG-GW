load gwc_output;

g0 = 0.067/pi;

%# get derivatives
dw = wss(2)-wss(1);
Nss = length(wss);
ddos = diff(dos)/dw;
ddos(Nss) = ddos(Nss-1);

dwc = wsc(2)-wsc(1);
Nsc = length(wsc);
ddosc = diff(dosc)/dwc;
ddosc(Nsc) = ddosc(Nsc-1);

%# smooth derivatives
width = 0.2; %# in meV
fss = exp(-wss.^2/width^2);
fss /= trapz(wss,fss);
fsc = exp(-wsc.^2/width^2);
fsc /= trapz(wsc,fsc);

sddos = dw*conv(fss,ddos);
sddos /= g0*27.21*1000
swss  = dw*[0:2*(Nss-1)]' + 2*wss(1);

sddosc = dwc*conv(fsc,ddosc);
swsc  = dwc*[0:2*(Nsc-1)]' + 2*wsc(1); 