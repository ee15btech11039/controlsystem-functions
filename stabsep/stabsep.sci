function [ga,gs,gi]=stabsep(g)
    //Decomposition of a system into stable and unstable part
    //
    //Calling Seqence
    //[ga,gs,gi]=stabsep(sys)
    //
    //Parameters
    //sys : a continuous or discrete time linear system.
    //ga  : a strictly proper antistable system.
    //gs  : a strictly proper stable system.
    //gi :  real matrix (or polynomial matrix for improper systems).
     
    //Description   
    //returns the stable-antistable decomposition of g.
    
    //Algorithm
    //Converting system into upper diagonal schur form i.e a=u'au,b=u'b,c=cu
    //such that a=|a11 a12 |  ,u : schur decomposition of a 
    //            |0   a22 |
    //further a=w^-1aw,b=w^-1b,c=cw
    //w=|I S | where S is solution of a11S-Sa22+a12=0
    //  |0 I |
    //decomposition of system : g=ga+gs
    
    //ref: http://www.incois.gov.in/proceedings/4.4%20backupofsksingsknagar.pdf
    
    //Examples
    //sys=ssrand(1,1,8);t=0.1;
    //sys=dscr(sys,t);
    //[ga,gs]=stabsep(sys);
 
    //Author
    //Ayush Kumar
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    [lhs,rhs]=argn(0),
    //error checking....................................................................................................................
    if rhs ~=1 then
        error(msprintf(gettext("%s: Wrong number of input arguments:%d argument expected.\n"),"stabsep",1))
    end,
    if and(typeof(g)<>['rational','state-space']) then
        error(msprintf(gettext("%s: Wrong type for input argument #%d: Linear state space or a transfer function expected.\n"),"stabsep",1))
    end,
    flag1=0;
    if (typeof(g)=='rational') then
        if(degree(g.num)>degree(g.den))
            error(msprintf(gettext("The %s command cannot be used for models with more zeroes than poles","stabsep")))
        end;
        g=tf2ss(g);
        flag1=1;
    end,
    if g.dt==[] then
        warning(msprintf(gettext("%s: Input argument %d is assumed continuous time.\n"),"stabsep",1));
        g.dt='c'
    end;
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //.....................................continuous system...............................................................//
    if(g.dt=='c') then
        [sla,sls,d]=dtsi(g); //stable-antistable decomposition for a 
        ga=sla;gs=sls;       //for a continuous time system
        if(flag1) then
                ga=clean(ss2tf(ga)),
                gs=clean(ss2tf(gs)),
        end
        return;
    end
    
    //......................................discrete system.................................................................//
    [a,b,c,d]=abcd(g),gi=d, //discrete system
    [n1,n2,t]=size(g),
    [a,u0]=balanc(a);b=u0\b;c=c*u0;
    [u,n]=schur(a,'d'),
    a=u'*a*u,  //first stage of transformation
    flag=0;
    if n==t then
        ga=0,
        flag='stable';
        gs=g,
    end,
    if n==0 then
        gs=0,
        ga=g,
        return;
    end,
    if(g.dt~='c') then
        if flag=='stable' then
            ga=0,gs=g;
            ga(7)=g(7);gs(7)=g(7);
            if(flag1) then
                ga=clean(ss2tf(ga)),
                gs=clean(ss2tf(gs));
            end
            return;
        
        else
    //      [ab,w,bs]=bdiag(a);
            a1=a(1:n,1:n),a4=a(n+1:t,n+1:t),x=a(1:n,n+1:t),
            z=sylv(a1,-a4,-x,'c'), //solving: a1z-za4+x=0
            w=[eye(n,n),z;0*ones(t-n,n),eye(t-n,t-n)], //second stage transformation
            wi=[eye(n,n),-z;0*ones(t-n,n),eye(t-n,t-n)],
            tr=u*w,tri=wi*u';
            bb=tri*b,b1=bb(1:n,:),b2=bb(n+1:t,:),
            cc=c*tr,c1=cc(:,1:n),c2=cc(:,n+1:t),
            ga=syslin('d',a4,b2,c2),ga(7)=g(7); //antistable part
            gs=syslin('d',a1,b1,c1);gs(7)=g(7); //stable part
            if(flag1) then  //returning system in rational form 
                ga=clean(ss2tf(ga)),
                gs=clean(ss2tf(gs));
            end
            return;
        end
        
    end;
endfunction

