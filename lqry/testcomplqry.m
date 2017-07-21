load('testlqqry.mat')
sys1=ss(a1,b1,c1,d1);%continuous
sys2=c2d(sys1,t1);%discrete 
[k1,x1]=lqry(sys1,q1,r1,n1')
[k2,x2]=lqry(sys1,q1,r1,n1')