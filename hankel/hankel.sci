function [nk,W]=hankel(g,tol)
    //Decomposition of a system into stable and unstable part
    //
    //Calling Seqence
    //[nk,W]=hankel(g,tol)
    //
    //Parameters
    //g   : a continuous or discrete time linear system.
    //tol : tolerance parameter default value is 1000*%eps)
    //nk  : squared hankel singular values 
    //W   : product of grammians .
     
    //Description   
    //returns the squared hankel singula values nk and W=P*Q where 
    //P=controllability grammian and Q=observability grammian.
    
    //Algorithm
    //decompose system to stable antistable part. g=ga+gs,
    //obtaining controllability and observebility grammians of gs.
    //nk:eigenvalues of product of grammians.
    
    //ref: http://web.iitd.ac.in/~janas/courses/material/eel879/sp_topics_01.pdf
    //http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=261486
    
    //Examples
    //sys=ssrand(1,1,8);t=0.1;
    //sys=dscr(sys,t);
    //nk=hankel(sys);
 
    //Author
    //Ayush Kumar
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    [lhs,rhs]=argn(0),
    //error checking .............................................................................................
    if rhs <1 || rhs >2  then
        error(msprintf(gettext("%s: Wrong number of input arguments:%d or %d argument expected.\n"),"hankel",1,2))
    end,
    if (typeof(g)=='rational') then
        if(degree(g.num)>degree(g.den))
            error(msprintf(gettext("The %s command cannot be used for models with more zeroes than poles","hankel")))
        end;
    end,

    if and(typeof(g)<>['rational','state-space']) then
        error(msprintf(gettext("%s: Wrong type for input argument #%d: Linear state space or a transfer function expected.\n"),"hankel",1))
    end
    
    if rhs==1 then 
    tol=100*%eps,
    end,
    
    if typeof(g)=='rational' then
        g=tf2ss(g); 
    end;
    lf=spec(g(2)),
    
    if min(abs(lf))<=tol then
        error(msprintf(gettext("%s: Wrong value for input argument #%d: Pure imaginary poles unexpected.\n"),"hankel",1)),
    end
    
    if g.dt==[] then
        warning(msprintf(gettext("%s: Input argument %d is assumed continuous time.\n"),"hankel",1));
        g.dt='c'
    end;
//................................continuous system...................................................................//
    if(g.dt=='c') then
        [sla,sls,d]=dtsi(g);
        lc=ctr_gram(sls),lo=obs_gram(sls),W=lc*lo;
        nk=gsort(real(spec(W)));
        return;
    end
//................................discrete system..................................................................//
    [a,b,c,d]=abcd(g),gi=d,  //decomposition of system
    [n1,n2,t]=size(g),
    [a,u0]=balanc(a);b=u0\b;c=c*u0;
    [u,n]=schur(a,'d'),
    a=u'*a*u,
    flag=0;
    if n==t then
        ga=0,
        flag='stable';
        gs=g,
    end,
    if n==0 then
        gs=0,
        ga=g,
    end,
        
    if(g.dt~='c') then
        if flag=='stable' then
            lc=ctr_gram(g),lo=obs_gram(g),W=lc*lo;
            nk=gsort(real(spec(W)));
            return;
        
        else
            a1=a(1:n,1:n),a4=a(n+1:t,n+1:t),x=a(1:n,n+1:t),
            z=sylv(a1,-a4,-x,'c'),
            w=[eye(n,n),z;0*ones(t-n,n),eye(t-n,t-n)],
            wi=[eye(n,n),-z;0*ones(t-n,n),eye(t-n,t-n)],
            tr=u*w,tri=wi*u';
            bb=tri*b,b1=bb(1:n,:),b2=bb(n+1:t,:),
            cc=c*tr,c1=cc(:,1:n),c2=cc(:,n+1:t),
            ga=syslin('d',a4,b2,c2),  //antistable part system
            gs=syslin('d',a1,b1,c1);  //stabel part system
            lc=ctr_gram(gs);
            lo=obs_gram(gs);
            W=lc*lo; //product of grammians
            nk=gsort(real(spec(W)));
            nk=clean(nk)
            return;
        end
        
    end;
endfunction

