function [K,X,e]=lqrd(A,B,Q,R,N,ts)
    //if typeof(P12)<>"state-space" then
      //  error(msprintf(gettext("%s: Wrong type for input argument #%d: Linear state space expected.\n"),"lqr",1))
    //end

    //if P12.dt == [] then
      //  error(msprintf(gettext("%s: Wrong value for input argument #%d: Time domain must be ''c'' or ''d''.\n"),"lqr",1))
    //end
    [lhs,rhs]=argn(0);
    [nx]=size(A,1)
    [nu]=size(B,2)
    //if rhs==5 then
      //  ts=N'
       // N=zeros(nu,nx)
    //end,
    
    if or(size(Q)<>nx) then
            error(msprintf(_("%s: Wrong size for input argument #%d: %d-by-%d matrix expected.\n"),"lqr",2,nx,nx))
    end

    if or(size(R)<>nu) then
            error(msprintf(_("%s: Wrong size for input argument #%d: %d-by-%d matrix expected.\n"),"lqr",3,nu,nu))
    end

    if norm(Q.'-Q,1)>100*%eps*norm(Q,1) then
            error(msprintf(_("%s: Wrong value for input argument #%d: Must be symmetric.\n"),"lqr",2))
    end

    if norm(R.'-R,1)>100*%eps*norm(R,1) then
            error(msprintf(_("%s: Wrong value for input argument #%d: Must be symmetric.\n"),"lqr",3))
    end
    if rhs~=6 then
        error(msprintf(_("%s: Wrong size for input argument #%d: %d-by-%d matrix expected.\n"),"lqr",4,nu,nx))
    end

    Q=(Q+Q')/2
    R=(R+R')/2
    vr = real(spec(R));
    vqnr = real(spec([Q N';N R]));
    if min(vr)<=0 then
        error(msprintf("%s Bad conditioning \n","lqrd"))
    elseif min(vqnr)<-100*%eps*max(0,max(vqnr)) then 
        warning(msprintf(gettext("%s Bad conditioning \n",'lqrd')))
    end
    
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
    QQ = (QQ+QQ')/2;        // Make sure QQ is symmetric
    Qd = QQ(1:nx,1:nx);
    Rd = QQ(nx+1:n,nx+1:n);
    Nd = QQ(1:nx,nx+1:n);
    Nd=Nd';
    ad = phi22(1:nx,1:nx);
    bd = phi22(1:nx,nx+1:n);
    
    [n]=size(B,1);
    I=eye(ad);Z=0*I;
    bigE=[I,Z,0*bd; ...
          Z,ad',0*bd; ...
          0*bd',-bd',0*bd'*bd];

    bigA=[ad,Z, bd; ...
         -Qd ,I, -Nd'; ...
          Nd, 0*bd', Rd];
            
   [w,ks]=schur(bigA,bigE,'d');
    if ks<>n then error(msprintf('lqrd: stable subspace too small!'));end
    ws=w(:,1:n);
    X12=ws(1:n,:);
    phi12=ws(n+1:2*n,:);
    u12=ws(2*n+1:2*n+nu,:);
    if rcond(X12)< 1.d-5 then warning('lqrd: bad conditionning!');end
    K=u12/X12;
    X=phi12/X12;
    K=-K
endfunction

