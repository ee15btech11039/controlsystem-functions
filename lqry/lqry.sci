function [K,X,e]=lqry(sys,Q,R,N)
    //Linear-quadratic compensator design with output weighting.
    //
    //Calling Seqence
    //output=lqry(sys,Q,R,N)
    //
    //Parameters
    //sys : state-space representation of a continuous or discrete time linear system.
    //Q : Real symmetric matrix,with same dimensions as sys.a
    //R : Full rank real symmetric matrix
    //N : real matrix, the default value is zeros(size(R,1),size(Q,2))
    //K : a real matrix, the optimal gain
    //X : a real symmetric matrix, the stabilizing solution of the Riccati equation
    //e : eigenvalues of sys.a-K*sys.b
     
    //Description   
    //[K,X,e] = lqry(sys,Q,R,N) 
    //Computes the linear optimal LQ full-state gain K for the linear dynamical system , 
    //the Riccati solution S, and the closed-loop eigenvalues e = eig(A-B*K).
    //          
    //  sys :     dx/dt=Ax+Bu
    // in continuous time or
    //          x+=Ax+Bu
    // in discrete time.
    //P:
    //   And instantaneous cost function X in l2-norm: 
    //      [y' u']BigQ[y u]' where BigQ=[Q N';N R]

    //Algorithm 
    //the function lqry is equivalent to lqr with weighting matrices:
    //  [QQ,NN';NN,RR]=[C',0;D',I]BigQ[C,D;0,I]
    
    //Examples
    //sys1=ssrand(3,1,3); //randomly generated sysytem
    //[a1,b1,c1,d1]=abcd(sys1);
    //[nx1,nu1]=size(sys1.b);
    //q1=rand(nx1,nx1);q1=(q1+q1')/2;
    //r1=rand(nu1,nu1);r1=(r1+r1')/2;
    //n1=rand(nu1,nx1);
    //[k1,x1]=lqry(sys1,q1,r1,n1)
    
    //Author
    //Ayush Kumar

    [lhs,rhs]=argn(0),
    //error checking ....................................................................
    if (rhs<3 || rhs>4) then
        error(msprintf(gettext("%s: Wrong number of input arguments:%d arguments sys,Q,R,N expected.\n"),"lqry",1))
    end
    if typeof(sys)<>"state-space" then
        error(msprintf(gettext("%s: Wrong type for input argument #%d: Linear state space expected.\n"),"lqry",1))
    end
    if(rhs==3)then
        N=zeros(size(R,1),size(Q,2)),
    end
    
    if sys.dt == [] then
        warning(msprintf(gettext("%s: Input argument %d is assumed continuous time.\n"),"lqry",1));
        g.dt='c'
    end
    A=sys.A
    B2=sys.B
    [nx,nu]=size(B2)
    if or(size(Q)<>nx) then
        error(msprintf(gettext("%s: Wrong size for input argument #%d: %d-by-%d matrix expected.\n"),"lqry",2,nx,nx))
    end
    if or(size(R)<>nu) then
        error(msprintf(gettext("%s: Wrong size for input argument #%d: %d-by-%d matrix expected.\n"),"lqry",3,nu,nu))
    end
    if norm(Q.'-Q,1)>100*%eps*norm(Q,1) then
        error(msprintf(gettext("%s: Wrong value for input argument #%d: Must be symmetric.\n"),"lqry",2))
    end
    if norm(R.'-R,1)>100*%eps*norm(R,1) then
        error(msprintf(gettext("%s: Wrong value for input argument #%d: Must be symmetric.\n"),"lqry",3))
    end,
    if argn(2)<4 then
        N=zeros(nu,nx);
    elseif or(size(N)<>[nu,nx]) then
        error(msprintf(gettext("%s: Wrong size for input argument #%d: %d-by-%d matrix expected.\n"),"lqry",4,nu,nx))
    end,
    [a,b,c,d]=abcd(sys)
    k2=[Q,N';N,R]
    if norm(k2.'-k2,1)>100*%eps*norm(k2,1) then
        warning(msprintf(gettext("%s:[Q,Nt;N,R]should be symmetric.\n"),"lqry"));
    end,
    row=size(k2)(1)
    col=size(k2)(2)
    lenr_Q=size(Q)(1)
    lenc_Q=size(Q)(2)
    lenr_c=size(c)(1)
    lenc_c=size(c)(2)
    zeromat=zeros(col-lenr_c,lenc_c)
    I=eye(col-lenr_c,size(d)(2))
    k1=[c',zeromat';d',I]
    k3=[c,d;zeromat,I]
    k=k1*k2*k3
    if norm(k.'-k,1)>100*%eps*norm(k,1) then
        warning(msprintf(gettext("%s: Bad conditionning.\n"),"lqry"));
        //warning(msprintf("%s:[C D;0 I]transx[Q,Ntran;N,R]x[C D;0 I]should be symmetric.\n"),"lqry"))
    end,
    q=k(1:lenr_Q,1:lenc_Q)
    n=k(1:lenr_Q,lenc_Q+1:size(k)(2))
    n=n'
    r=k(lenr_Q+1:size(k)(1),lenc_Q+1:size(k)(2))
    [K,X]=lqr(sys,q,r,n)
    e=spec(a+b*K)
    K=-K
    return;
endfunction
