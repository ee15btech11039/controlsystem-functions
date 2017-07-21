//test case
sys=ssrand(1,1,4)
[a,b,c,d]=abcd(sys);t=0.1;
sysdscr=dscr(sys,t) //discrete unstable system
savematfile('testhankel.mat','a','b','c','d','t');
n1=hankel(sys);
disp("n1 : continuos system")
disp(n1^0.5)
n2=hankel(sysdscr);
disp("n2 : discrete system")
disp(n2^0.5)
