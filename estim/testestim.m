load('testestim.mat')
sys1=ss(a1,b1,c1,d1);
sys2=ss(a2,b2,c2,d2);
rsys1=estim(sys1,l1);
rsys2=estim(sys2,l2,sensors,known);
rsys1.a,rsys1.b,rsys1.c,rsys1.d
rsys2.a,rsys2.b,rsys2.c,rsys2.d