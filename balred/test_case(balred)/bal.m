sys1=rss(10)
a1=sys1.a
b1=sys1.b
c1=sys1.c
d1=sys1.d


systest=rss(4,2) %stable system
aa=[-0.5566,1.0000,-1.2472,-0.8536;-0.0197,-0.5566,0.8283,0.5669;0,0,1.7550,1.4611;0,0,0,0.8605]
bb=[0,0,0,2.9837]'
cc=[3.0574,0,1.0898,0.7459]
dd=0
sys3_=ss(aa,bb,cc,dd) %unstable system
sys3=sys3_+systest   %unstable system
[a,b,c,d]=ssdata(sys3)
save('sys10.mat','a1','b1','c1','d1','a','b','c','d')
