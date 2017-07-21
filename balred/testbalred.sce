//test case

sys=ssrand(1,1,10); //a random 10th order unstable system
[a,b,c,d]=abcd(sys);
savematfile("checkbalred.mat",'a','b','c','d')
k=balred(sys,6);
ff=0.01:0.01:10;
bode(sys,ff);
title("original system")
scf();
bode(k,ff);
title("reduced system")
