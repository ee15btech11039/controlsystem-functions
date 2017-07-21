
function[u]=filtbutter(data,n,band)
    //nth order butterworth  filter
    //
    //Calling Seqence
    //[u]=filtbutter(data,n,band)
    //
    //Parameters
    //data : input signal
    //n : order of the butterworth filter
    //band: band of the filter
     
    //Description   
    //[u]=filtbutter(data,n,band) returns the filtered input signal using an nth order 
    //butterworth filter
    
    //Examples
    // data=rand(50,2)
    //u=filtbutter(data,8,[0.1 0.9])
    
    //Author
    //Ayush Kumar
    
    //////////////////////////////////////////////////////////////////////////////////////////////////
    
    [lhs,rhs]=argn(0),
    [nx,nu]=size(data);
    //error checking....................................................
    if rhs~=3 then
        error(msprintf(gettext("%s wrong number of input argument %d arguments expected","filt",2)))
    end,
    if nx<=3*n then
        error(msprintf(gettext("%s number of samples in data must be atlest 3 times order of filter","filt")))
    end,
    if and(size(band)==[1,1]) then
        band=[0,band];
    end
    if and(size(band)~=[1 2]) then
        error(msprintf(gettext("%s :wrong type of input argument %d","filt",3)))
    end
    if band(2)<band(1) then
        error(msprintf(gettext("%s :band first elemnt should be less than second element"),"filt")),
    end
    if typeof(n)~="constant" || or(size(n)~=[1,1]) then
        error(msprintf(gettext("%s :wrong type of input argument %d"),"filt",2)),
    end
    band=band/2  
    if ~and(band==[0 0.5]) then
        if(band(1)==0) then
            [hz]=iir(n,'lp','butt',[band(2) band(1)],[0 0]);  //nth order analog butterwoth filter
            num=hz(2);
            den=hz(3);
            for i=1:1:nu
              y(:,i)=filter(num,den,data(:,i));
            end
        elseif(band(2)==0.5) then
            [hz]=iir(n,'hp','butt',[band(1) band(2)],[0 0]);  
            num=hz(2);
            den=hz(3);
            for i=1:1:nu
                y(:,i)=filter(num,den,data(:,i));
            end
        else
            [hz]=iir(n,'bp','butt',band,[0 0]);  
            num=hz(2);
            den=hz(3);
            for i=1:1:nu   //digital filter
                y(:,i)=filter(num,den,data(:,i));
            end
        end
            
        u=y
    
    else
        u=data;
    end
        
endfunction



