//example1
sys1=ssrand(1,1,2);
[a1,b1,c1,d1]=abcd(sys1);
k1=[1 3];
l1=[2 5]';
[ac1,bc1,cc1,dc1]=reg(sys1,k1,l1) //returns state spaces of the regulator


//example2
sys2=ssrand(3,3,4);
[a2,b2,c2,d2]=abcd(sys2);
k2=[1 2 4 3;3 4 1 2] ; 
l2=[1 4 5 4]';
controls=[1 2];
sensors=[3]
known=[3]
savematfile('testreg.mat','a1','b1','c1','d1','k1','l1','a2','b2','c2','d2','k2','l2','controls','sensors','known')
[ac2,bc2,cc2,dc2]=reg(sys2,k2,l2,sensors,known,controls)

