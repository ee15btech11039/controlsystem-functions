function[rsys]=balred(sys,n)
    //a reduced-order approximation of the input LTI system
    //
    //Calling Seqence
    //output=balred(sys,n)
    //
    //Parameters
    //sys : state-space or rational model of a continuous or discrete time linear system.
    //n : the desired reduced order of the sytem,can be a number or vector or matrix
    //rsys : reduced order approximated system, according to the input n
     
    //Description   
    //[rsys] = balred(sys,n) 
    //computes the reduced order approximation of a given lti systetm,when sys has 
    //unstable poles ,it is first decomposed into its stable and unstable part
    //            sys=ga+gs              ,and only 
    //the stable part is approximated,after approximation the final system system is 
    //returned as the sum of reduced stable and unstable part.
    //   
    //Algorithms
    //blanced truncation algorithm is used to reduce the order of the system.
    //
    //ref: http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=261486&tag=1
    //
    //Examples
    //sys=ssrand(1,1,10);
    //sysred=balred(sys,6);
    //ff=0.01:0.01:100;
    //bode(sys,ff);
    //scf(1);
    //bode(sysred,ff);
    //
    //Author
    //Ayush Kumar
    
    
    [lhs,rhs]=argn(0),
    //................................................................................
    //error checking
    if rhs ~=2 then
        error(msprintf(gettext("%s: Wrong number of input arguments:  %d arguments expected.\n"),"balred",2))
    end,
    if and(typeof(sys)<>['rational','state-space']) then
        error(msprintf(gettext("%s: Wrong type for input argument #%d: Linear state space or a transfer function expected.\n"),"balred",1))
    end 
    if typeof(n)~='constant' then
        error(msprintf(gettext("%s: Wrong type for input argument #%d: positive integer or vectors expected\n"),"balred",2))
    end ,
    if ~and(n>=0) then
        error(msprintf(gettext("%s: Wrong type for input argument #%d: positive integer or vectors expected\n"),"balred",2))
    end , 
    if ~and(n-ceil(n)==0) then
        error(msprintf(gettext("%s: Wrong type for input argument #%d: positive integer or vectors expected\n"),"balred",2))
    end , 
     
    if sys.dt==[] then
        warning(msprintf(gettext("%s: Input argument %d is assumed continuous time.\n"),"balred",1));
        sys.dt='c'
    end;
    
    if sys.dt=='c' then
        K=0
    else
        K=1
    end
    
    
    //converting  rational form to state space form if required.........................
    flag=0;
    if typeof(sys)=='rational' then
        if(degree(sys.num)>degree(sys.den))
            error(msprintf(gettext("The %s command cannot be used for models with more zeroes than poles","balred")))
        end;
        flag='rational',
        sys=tf2ss(sys),
    end
    
    [nx,nu]=size(n),
    
    len_n=nx*nu, //size of the rsys array in which we will store the data
    //if len_n>1 then
     //   rsys=cell(len_n,1),
    //end
    
    //.................................................................................
    order=size(sys(2)) 
    order=order(1) // order of the original system
    for i=1:nx
        for j=1:nu
            if n(i,j)>order then
                error(msprintf(gettext("%s:wrong input argument:%d argument should be less than or equal to order of the system","balred",2)))
            end
        end
    end
    
    //
    //order reduction can be done for the stable part only unstable part of the system will remain unaffeted.
    //Decomposition of the system to stable and unstable part
    
    [ga,gs]=stabsep(sys);
    flag2=0;
    if(ga==0)then
        flag2=1,
    end,
    if(gs==0)then  //if system is completely unstable it can't be reduced
        warning(msprintf(gettext("%s:as stable part of the system is constant system cannot be reduced,returning the original system \n"),"balred"));
        if(len_n==1)then
            rsys=ga;
            return;
        else
            for(i=1:len_n)
                rsys{i,1}=ga;
            end
            return;
        end,
    else
        gs=balreal(gs);
    end,
     //here we are using the balaned realization truncation algorithm
    [a,b,c,d]=abcd(gs); 
    len=size(a,1) //order of the stable part of the system
    count=1
    for i=1:nx
        for j=1:nu
            
            orderred=order-n(i,j)
            if orderred>len then//warning message type and reduce the stable part to 0th order
                warning(msprintf(gettext("%s:the approximation order should be greater than or equal to number of unstable poles\n"),"balred"));
                orderred=len
            end
            nn=len-orderred //order to which we have to reduce the stable part of the system
    
    //balanced truncation algorithm implementation.......
            if nn<len & nn~=0 then
                a11=a(1:nn,1:nn)
                a12=a(1:nn,nn+1:len)
                a21=a(nn+1:len,1:nn)
                a22=a(nn+1:len,nn+1:len)
                b1=b(1:nn,:)
                b2=b(nn+1:len,:)
                c1=c(:,1:nn)
                c2=c(:,nn+1:len)
                I=eye(size(a22)(1),size(a22)(2))
                A=a11+a12*inv(K*I-a22)*a21
                B=b1+a12*inv(K*I-a22)*b2
                C=c1+c2*inv(K*I-a22)*a21
                D=d+c2*inv(K*I-a22)*b2
            end
            if nn==len & nn~=0 then
                A=a,B=b,C=c,D=d
        
            end
            if nn==0 then
                a11=0,a12=0,a21=0,a22=a
                b1=0,b2=b
                c1=0,c2=c
                I=eye(size(a22)(1),size(a22)(2))
                A=a11+a12*inv(K*I-a22)*a21
                B=b1+a12*inv(K*I-a22)*b2
                C=c1+c2*inv(K*I-a22)*a21
                D=d+c2*inv(K*I-a22)*b2
        
            end
    
            
            if sys.dt=='c' then
                gs=syslin('c',A,B,C,D)
                if (flag2) then
                    rsysl=gs,
                else
                    [lenr_ga,lenc_ga]=size(ga.a)
                    rsysl=ga+gs,
                    f_order=nn+lenr_ga,
                    rsysl.a=rsysl.a(1:f_order,1:f_order)
                    rsysl.b=rsysl.b(1:f_order,:)
                    rsysl.c=rsysl.c(:,1:f_order)
                    rsysl(6)=rsysl(6)(1:f_order,:)
                end,
                if flag=='rational' then
                    rsysl=ss2tf(rsysl)
                    rsysl=clean(rsysl)
                end
                if len_n==1 then
                    rsys=rsysl,
                    return;
                else
                    rsys{count,1}=rsysl
                    count=count+1,
                end,
            else
                gs=syslin('d',A,B,C,D)
                gs(7)=sys(7)
                if (flag2) then
                    rsysl=gs,
                else
                    [lenr_ga,lenc_ga]=size(ga.a)
                    rsysl=ga+gs,
                    f_order=n+lenr_ga,
                    rsysl.a=rsysl.a(1:f_order,1:f_order)
                    rsysl.b=rsysl.b(1:f_order,:)
                    rsysl.c=rsysl.c(:,1:f_order)
                    rsysl(6)=rsysl(6)(1:f_order,:)
                end,
                if flag=='rational' then
                    rsysl=ss2tf(rsysl)
                    rsysl=clean(rsysl)
                end
                if len_n==1 then
                    rsys=rsysl,
                    return;
                else
                    rsys{count,1}=rsysl
                    count=count+1,
                end,
            end
        end
    end
endfunction
 
    
    
    
