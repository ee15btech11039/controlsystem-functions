
u1=idinput(100) //100 length single channel signal
u2=idinput([100 2 2]) //2 channel signal of length 200 and period 100
u3=idinput([100 2 2],'rbs',[0.1 0.65])//signal passed through band [0.1 0.65]
u4=idinput([100 2 2],'rbs',[0.1 0.65],[-1 3])// changing the levels
