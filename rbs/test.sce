//filter result 
//here in test case we are comparing the outputs of idfilt and filt function defined 
//here for a given input ,if the output of these filters will be same for given same 
//input the orresponding result of rbs function will also be same as rest all part of
//rbs function just adjusts the outputs from the filter  with their corresponding sign 
//and levels. 

loadmatfile('testfile2.mat')

function[u]=filt(data,band1,nu)
    [lhs,rhs]=argn(0)
    band=band1/2  //dividing the band by 2 as we are taking data from matlab and here band is defined in between [0 0.5].while in matlab it is defined in bw [0 1]
    if ~and(band==[0 0.5]) then
        if(band(1)==0) then
            [hz]=iir(8,'lp','butt',[band(2) band(1)],[0 0]);  //8th order butterwoth filter
            num=hz(2);
            den=hz(3);
            for i=1:1:nu
              y(:,i)=filter(num,den,data(:,i));
            end
        elseif(band(2)==0.5) then
            [hz]=iir(8,'hp','butt',[band(1) band(2)],[0 0]);  //8th order butterwoth filter
            num=hz(2);
            den=hz(3);
            for i=1:1:nu
                y(:,i)=filter(num,den,data(:,i));
            end
        else
            [hz]=iir(8,'bp','butt',band,[0 0]);  //8th order butterwoth filter
            num=hz(2);
            den=hz(3);
            for i=1:1:nu
                y(:,i)=filter(num,den,data(:,i));
            end
        end
            
        u=y
    
    else
        u=data;
    end
        
endfunction

testu1=filt(u1,band11,nu1);
testu2=filt(u2,band21,nu2);
testu3=filt(u3,band31,nu3);
testu4=filt(u4,band41,nu4);
testu5=filt(u5,band51,nu5);

