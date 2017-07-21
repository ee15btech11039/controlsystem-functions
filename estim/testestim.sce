//example
sys1=ssrand(1,1,2)
[a1,b1,c1,d1]=abcd(sys1);
l1=[1 2]'
savematfile('testestim.mat','a','b','c','d')
[ae1,be1,ce1,de1]=estim(sys1,l1)


//Example
sys2=ssrand(2,2,2);
[a2,b2,c2,d2]=abcd(sys2);
l2=[1 2;3 4]
sensors=[1 2]
known=[2]
[ae2,be2,ce2,de2]=estim(sys2,l2,sensors,known)

savematfile('testestim.mat','a1','b1','c1','d1','l1','a2','b2','c2','d2','l2','sensors','known')
