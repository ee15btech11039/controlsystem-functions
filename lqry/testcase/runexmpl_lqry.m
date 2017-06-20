load('testlqry4.mat')
sys2=c2d(sys1,0.01) %discrete system
[k1,x1]=lqry(sys1,q1,r1,n1)
[k2,x2]=lqry(sys2,q1,r1,n1)


%sys2=c2d(sys2,0.5)
%sys3=c2d(sys3,0.01)


%[k3,x3]=lqry(sys3,q3,r3,n3)
 