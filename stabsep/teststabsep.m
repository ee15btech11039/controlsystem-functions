load('checkstabsep.mat');
sys=ss(a,b,c,d,t);
[gs,ga]=stabsep(sys);
opts = bodeoptions('cstprefs');
opts.FreqUnits = 'Hz';
figure(1);
bodeplot(ga,opts);
title('antistable part')
figure(2);
bode(gs,opts);
title('stable part');
