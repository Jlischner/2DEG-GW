%# input file
n_cm = 10*10^10; %# density in 1/cm^2
d = 230/0.5291; %# size of quantum well
d_el = 245/0.5291; %# distance from center of well to electrode
nus = [-10: 0.1 : 10]'; %# freq. grid for sigma in units of EF
#kys = [0:0.01:1]'; %# kpoint in units of kF
kys = [0:1]';
global m  = 0.067; %# Ashoori: 0.067
global eps_back = 10.9;
Lambda = 20; %$ kc/kF
etaEF = 0.05; %# broadening for correlation self energy