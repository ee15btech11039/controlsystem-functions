data=rand(100,3);
n=5;band=[0.2 0.8];
u=filtbutter(data,n,band)
savematfile("testbutter.mat",'data','n','band')


//here for a given input ,if the output of these filters will be same for given same 
//input the corresponding result of rbs function will also be same as rest all part of
//rbs function just adjusts the outputs from the filter  with their corresponding sign 
//and levels. 
