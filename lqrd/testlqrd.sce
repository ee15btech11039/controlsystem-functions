sys1=ssrand(3,1,3); //randomly generated sysytem
a1=sys1.A;b1=sys1.B;

[nx1,nu1]=size(sys1.b);
q1=rand(nx1,nx1);q1=(q1+q1')/2;
r1=rand(nu1,nu1);r1=(r1+r1')/2;
n1=rand(nu1,nx1);
t1=0.1
[k11,x1]=lqrd(a1,b1,q1,r1,n1,t1) //%continuos

savematfile('testlqrd.mat','a1','b1','q1','r1','n1','t1')

