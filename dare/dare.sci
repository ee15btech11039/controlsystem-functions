
function [X,L,G]=dare(a,b,q,varargin)
    //solving general discrete algebric riccati equation
    //
    //Calling Seqence
    //[X,L,G]=dare(a,b,q)
    //[X,L,G]=dare(a,b,q,r)
    //[X,L,G]=dare(a,b,q,r,s)
    //
    //Parameters
    //a : real matrix (n-by-n).
    //b : real symmetric matrix (n-by-m).
    //q : real symmetric matrix (n-by-n).
    //r : real matrix (m-by-m).default(eye(m,m))
    //s :real matrix (n-by-m)
    //X : unique stabilizing solution of the discrete-time algebric Riccati equation.
    //L : closed-loop poles.
    //G : corresponding gain matrix.
     
    //Description   
    //[X]=dare(A,B,Q,R) computes the unique stabilizing solution X of the 
    //discrete-time algebraic Riccati equation
    // A'XA-X-A'XB(B'XB+R)^-1(B'XA)+Q=0
    //[X] =dare(A,B,Q,R,S) solves the more general discrete-time algebraic 
    //Riccati equation
    //A'XA-X-(A'XB+S)(B'XB+R)^(-1)(B'XA+S')+Q=0
    //
    //Algorithm
    //The general solution of the riccati equation is obtained by schur 
    //factorisation of the matrix pencils associated with these Riccati 
    //equations :       
    //z[I BR^(-1)B' ; 0  (A-BR^(-1)S')'] - [A-BR^(-1)S' 0 ; SR^(-1)S'-Q].
    //
    //ref: http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=1457351&tag=1
    //   
    //Example
    //a=rand(4,4)
    //b=rand(4,2)
    //q=rand(4,4);q=(q+q')/2;
    //r=rand(2,2);r=(r+r')/2;
    //s=rand(4,2);
    //x=dare(a,b,q,r,s);
    //disp(x)
    //
    //Author
    //Ayush Kumar
    
    
    [lhs,rhs]=argn(0),
    //error checking for dimensions
    if rhs<3 || rhs>5  then
        error(msprintf(gettext("%s:wrong number of input arguments","dare")))
    end,
    
    [nx,nu]=size(a);
    [nxb,nub]=size(b)
    if nx~=nu then
        error(msprintf(gettext("%s: Wrong size for input argument #%d: a square matrix expected.\n"),"dare",1));
    end,
    if or(size(q)~=nx) then
            error(msprintf(gettext("%s: Wrong size for input argument #%d: %d-by-%d matrix expected.\n"),"dare",3,nx,nx))
    end,
    if (size(b,1)~=nx) then
            error(msprintf(gettext("%s: Wrong size for input argument #%d: its number of rows should be same as %s.\n"),"dare",2,"a"))
    end,
    if rhs==3 then
        r=eye(nub,nub)  //a identity square matrix having same columns as of b  
    end
    if rhs==4 then
        r=varargin(1)
    end,
    if or(size(r)~=nub)then
        error(msprintf(gettext("%s: Wrong size for input argument #%d: %d-by-%d matrix expected.\n"),"dare",3,nx,nx))
    end,
    if norm(r.'-r,1)>100*%eps*norm(r,1) then
            error(msprintf(_("%s: Wrong value for input argument #%d: Must be symmetric.\n"),"dare",4))
    end
    s=0*b;
    if rhs==3||rhs==4 then
        //dare(a,b,q,r)
        bb=b/r*b'
        X=riccati(a,bb,q,'d') 
        G=inv(b'*X*b+r)*(b'*X*a+s')
        L=spec(a-b*G)
        return;
    end,
    
    if rhs==5 then
        r=varargin(1);
        s=varargin(2);
    end,
    
    if size(s,1)~=nxb || size(s,2)~=nub then
        error(msprintf(gettext("%s: Wrong size for input argument #%d: %d-by-%d matrix expected.\n"),"dare",5,nxb,nub))
    end,
    
    [n,nru]=size(b);
    
    I=eye(a);Z=0*I;Ri=inv(r)
    
    aa=[I b*Ri*b';0*ones(n,n) (a-b*Ri*s')'];//pencils associated with the riccati equation
    bb=[a-b*Ri*s'  0*ones(n,n);s*Ri*s'-q eye(nx,nx)'];
    
    [bs,as,ss,n1]=schur(bb,aa,"d");//schur factorisation of the pencils
        if n1<>n then
            error(msprintf(gettext("%s: Wrong dimension (%d) of stable subspace: %d expected.\n"),"dare",n1, n))
        end
        ss=ss(:,1:n1);
    
    x1=ss(n+1:2*n,:),x2=ss(1:n,:),
    X=x1/x2,
    G=inv(b'*X*b+r)*(b'*X*a+s'),
    L=spec(a-b*G),
    
endfunction

