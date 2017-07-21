load('testhankel.mat');
sys=ss(a,b,c,d);
n1=hsvd(sys)
sysdscr=c2d(sys,t);
n2=hsvd(sysdscr)
