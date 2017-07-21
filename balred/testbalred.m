load('checkbalred.mat');
sys=ss(a,b,c,d);
kmat=balred(sys,6);
opts = bodeoptions('cstprefs');
opts.FreqUnits = 'Hz';
bodeplot(kmat,opts);
title('reduced system')