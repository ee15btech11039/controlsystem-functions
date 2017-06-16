loadmatfile('sys10.mat')
sys1=syslin('c',a1,b1,c1,d1)  //continous system 10th order
sys2=dscr(sys1,0.1)  // discrete system
sys3=syslin('c',a,b,c,d)  //an 8th order unstable syatem

sys21_sci=balred(sys1,6)
sys22_sci=balred(sys2,4)
sys23_sci=balred(sys3,4)

[a21,b21,c21,d21]=abcd(sys21_sci)
t21=sys21_sci(7)

[a22,b22,c22,d22]=abcd(sys22_sci)
t22=sys22_sci(7)

[a23,b23,c23,d23]=abcd(sys23_sci)


savematfile('exbalred1.mat','a21','b21','c21','d21','t21','a22','b22','c22','d22','t22','a23','b23','c23','d23')
//sys21=clean(ss2tf(sys21))
//sys22=clean(ss2tf(sys22))
//
//disp("4th order reduced model",sys22)
