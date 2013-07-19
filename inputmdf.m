load("gwc_input.dat");
ws = dat(:,1);
ReSig = dat(:,2);
ImSig = dat(:,3);
vxc = 0;
eF = 0;
ts   = [-50000: 10 : 50000]';
wsN = [-0.03 :0.00005: 0.005]';
useEqp = 1; %# use Eqp in cumulant from Dyson eq.
zeroalpha = 0; %# set asymmetry to zero
fitrange = 5; %# deal with w=0 in C(w)
qp_thresh = 800; %# helps finding the qp-peak if the satellite peak is higher