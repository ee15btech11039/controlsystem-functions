sys1=rss(5,5)
[a1,b1,c1,d1]=ssdata(sys1);
[nx1,nu1]=size(sys1.b);
q1=randn(nx1,nx1);q1=q1'*q1; %q and r matrices should be positive semi definite
r1=randn(nu1,nu1);r1=r1'*r1;
n1=randn(nx1,nu1);
save('testlqry4.mat','a1','b1','c1','d1','q1','r1','n1')

