
loadmatfile('testlqry4.mat')
sys1=syslin('c',a1,b1,c1,d1)
[k1,x1]=lqry(sys1,q1,r1,n1')
sys2=dscr(sys1,0.01)
[k2,x2]=lqry(sys2,q1,r1,n1')


disp("k1")
disp(k1)
disp("x1")
disp(x1)

disp("k2")
disp(k2)
disp("x2")
disp(x2)

 
//[k2,x2]=lqry(sys2,q2,r2,n2')
