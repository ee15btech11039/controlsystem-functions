load('sys10.mat')%system loaded from here(sys1,sys3)
load('exbalred1.mat')%scilab results stored here
sys1=ss(a1,b1,c1,d1)  %example1
sys2=c2d(sys1,0.1)    %example2 discrete system
sys3=ss(a,b,c,d)

sys21=balred(sys1,6)%matlab results
sys21_sci=ss(a21,b21,c21,d21)%scilab result continous system
sys22=balred(sys2,4)
sys22_sci=ss(a22,b22,c22,d22,t22)
sys23=balred(sys3,4)
sys23_sci=ss(a23,b23,c23,d23)


figure(1);
bode(sys21)
legend('matlab result')%result comparision
title('redn of 10 order system to 6th order')
hold on
bode(sys21_sci,'r--')

figure(2);
bode(sys22)
legend('matlab result')
hold on
bode(sys22_sci,'r--')
title('redn of 10 order discrete system to 4th order')
figure(3);
bode(sys23)
legend('matlab result')
hold on
bode(sys23_sci,'r--')
title('redn of 8 order unstable system to 4th order')


