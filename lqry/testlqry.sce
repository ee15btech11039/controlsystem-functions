
sys1=ssrand(3,1,3); //randomly generated sysytem
[a1,b1,c1,d1]=abcd(sys1);
[nx1,nu1]=size(sys1.b);
q1=rand(nx1,nx1);q1=(q1+q1')/2;
r1=rand(nu1,nu1);r1=(r1+r1')/2;
n1=rand(nu1,nx1);
[k1,x1]=lqry(sys1,q1,r1,n1) //%continuos
t1=0.1
sys2=dscr(sys1,t1);
[k2,x2]=lqry(sys2,q1,r1,n1) //discrete
savematfile('testlqqry.mat','a1','b1','c1','d1','q1','r1','n1','t1')
