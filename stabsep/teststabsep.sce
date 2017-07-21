//test case 
sys=ssrand(1,1,10);t=0.1;
sys=dscr(sys,t);
[ga,gs]=stabsep(sys);
[a,b,c,d]=abcd(sys);
savematfile('checkstabsep.mat','a','b','c','d','t');
scf();
bode(ga);
title("antistabel part");
scf();
bode(gs);
title("stable part");
