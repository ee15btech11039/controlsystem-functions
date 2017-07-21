//test example
a=rand(4,4)
b=rand(4,2)
q=rand(4,4);q=(q+q')/2;
r=rand(2,2);r=(r+r')/2;
s=rand(4,2);
savematfile('checkdare.mat','a','b','q','r','s')
[x,l,g]=dare(a,b,q,r,s);

