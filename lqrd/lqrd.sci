function [Kd,Xd]=lqrd(A,B,Q,R,N,ts)
    //discrete linear-quadratic (LQ) regulator for continuous plant
    //
    //Calling Seqence
    //[Kd,Xd]=lqrd(A,B,Q,R,ts)
    //[Kd,Xd]=lqrd(A,B,Q,R,N,ts)
    //
    //Parameters
    //A,B: continuous plant state spaces
    //Q : real symmetric matrix,with same dimensions as A
    //R : real symmetric matrix,with same dimension as B
    //N : real matrix,the default value is zeros(size(B,2),size(B,1))
    //ts :sampling time of the discrete regulator
    
    //Description   
    //[K,X,e] = lqrd(sys,Q,R,N,ts) 
    // designs a discrete full-state-feedback regulator that has response 
    //characteristics similar to a continuous state-feedback regulator designed 
    //using lqr.
    //          

    //Algorithm 
    //The equivalent discrete gain matrix Kd is determined by discretizing the 
    //continuous plant and weighting matrices(Qd,Rd,Nd) using the sample time Ts and the 
    //zero-order hold approximation.
    //Xd is calculated by solving
    //A'XdA−Xd−(A'XdB+Nd)(B'XdB+Rd)^-1(B'XdA+Nd')+Qd=0 from which Kd is obtained as
    // Kd=(B'XdB+R)^-1(B'XdA+Nd')
    
    //Examples
    //sys1=ssrand(3,1,3); //randomly generated sysytem
    //a1=sys1.A;b1=sys1.B;
    //[nx1,nu1]=size(sys1.b);
    //q1=rand(nx1,nx1);q1=(q1+q1')/2;
    //r1=rand(nu1,nu1);r1=(r1+r1')/2;
    //n1=rand(nu1,nx1);
    //t1=0.1
    //[kd,xd]=lqrd(a1,b1,q1,r1,n1,t1) 
    
    //ref :http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=1101743
    
    //Author
    //Ayush Kumar
    
    [lhs,rhs]=argn(0);
    [nx]=size(A,1)
    [nu]=size(B,2)
    
     
    if and(rhs<>[5,6]) then
        error(msprintf(_("%s: wrong number of input arguments %d or %d arguments expected.\n"),"lqrd",5,6))
    end
    if rhs==5 then
        ts=N
        N=zeros(nu,nx) //default value of N
    end,
    if or(size(Q)<>nx) then
            error(msprintf(_("%s: Wrong size for input argument #%d: %d-by-%d matrix expected.\n"),"lqrd",3,nx,nx))
    end,

    if or(size(R)<>nu) then
            error(msprintf(_("%s: Wrong size for input argument #%d: %d-by-%d matrix expected.\n"),"lqrd",4,nu,nu))
    end,
    if or(size(N)<>[nu,nx]) then
            error(msprintf(_("%s: Wrong size for input argument #%d: %d-by-%d matrix expected.\n"),"lqrd",5,nu,nx))
    end,

    if norm(Q.'-Q,1)>100*%eps*norm(Q,1) then
            warning(msprintf(_("%s: Q matrixmust be symmetric.\n"),"lqrd"))
    end

    if norm(R.'-R,1)>100*%eps*norm(R,1) then
            warning(msprintf(_("%s: R matrix must be symmetric.\n"),"lqrd"))
    end
   
    //making Q and R matrices symmetric
    Q=(Q+Q')/2
    R=(R+R')/2
    //cheking positivity
    vr = real(spec(R));
    vqnr = real(spec([Q N';N R]));
    if min(vr)<=0 then
        error(msprintf(gettext("%s Bad conditioning \n","lqrd")))
    elseif min(vqnr)<-100*%eps*max(0,max(vqnr)) then 
        warning(msprintf(gettext("%s Bad conditioning \n","lqrd")))
    end
    //discrete equivalent of continuous cost function
    n=nx+nu
    Za = zeros(nx,nx); Zb = zeros(nx,nu); Zu = zeros(nu,nu);
    M = [ -A' Zb   Q  N';
          -B' Zu   N  R;
           Za  Zb   A   B;
           Zb' Zu  Zb' Zu];
    phi = expm(M*ts);
    phi12 = phi(1:n,n+1:2*n);
    phi22 = phi(n+1:2*n,n+1:2*n);
    QQ = phi22'*phi12;
    QQ = (QQ+QQ')/2;    //making QQ symmetric
    Qd = QQ(1:nx,1:nx);
    Rd = QQ(nx+1:n,nx+1:n);
    Nd = QQ(1:nx,nx+1:n);
    Nd=Nd';
    ad = phi22(1:nx,1:nx); //discretized sytem matrices
    bd = phi22(1:nx,nx+1:n);
    //solving discrete time algebric riccati equation.
    [n]=size(B,1);
    I=eye(ad);Z=0*I;
    bigE=[I,Z,0*bd; ...
          Z,ad',0*bd; ...
          0*bd',-bd',0*bd'*bd];

    bigA=[ad,Z, bd; ...
         -Qd ,I, -Nd'; ...
          Nd, 0*bd', Rd];
            
   [w,ks]=schur(bigA,bigE,'d');
    if ks<>n then error(msprintf(gettext("lqrd: stable subspace too small!")));end
    ws=w(:,1:n);
    X12=ws(1:n,:);
    phi12=ws(n+1:2*n,:);
    u12=ws(2*n+1:2*n+nu,:);
    if rcond(X12)< 1.d-5 then warning(msprintf(gettext("lqrd: bad conditioning!")));end
    Kd=u12/X12;
    Xd=phi12/X12;
    Kd=-Kd
endfunction
