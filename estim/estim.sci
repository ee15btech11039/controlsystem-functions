function [ae,be,ce,de,ts]=estim(sys,l,sensors,known)
    //Return state spaces of a state estimator for a given estimator gain.
    //
    //Calling Sequence
    //[ae,be,ce,de]=estim(sys,l)
    //[ae,be,ce,de]=estim(sys,l,sensors,known)
    //
    //Parameters
    //sys:lti model
    //l  :state feedback matrix.
    //sensors : Indices of measured output signals y,if omitted all outputs are measured.
    //known : Indices of known input signals u (deterministic) to sys. All other inputs to sys
    //are assumed stochastic.Default known=[]
    //[ae,be,ce,de]:state spaces of the estimator
    
    //Description   
    //estim(sys,l)produces a state estimator est given the plant state-space model sys and the 
    //estimator gain L.All inputs of sys are assumed stochastic and all outputs are measured
    // The estimator est state-spaces are returned.
    //For a continuous system  dx/dt=Ax+Bw,y=Cx+Dw
    //estim uses the following equation to generate plant output
    //           _                      _
    //  estimate y and a state estimate x.
    //  _     _      _
    // dx/dt=Ax+L(y-Cx)
    //   _
    // | y |  |C| _
    // | _ |= | | x     
    // | x |  |I|
    //
    //
    //Author
    //Ayush Kumar
    
    [lhs,rhs]=argn(0);
    [a,b,c,d]=abcd(sys);
    [nx]=size(a,1);
    [ny,nu]=size(d);
    //error checking
    if rhs<2 || rhs>4 then
        error(msprintf(gettext("%s : wrong number of input arguments "),"estim"))
    end,
    if typeof(sys)~="state-space" then
        error(msprintf(gettext("%s : state-space model expected "),"estim"))
    end
    if or(size(l)<>[nx,ny]) then
        error(msprintf(gettext("%s :l should have as many rows as states and as many columns as measured output "),"estim"))
    end
    
    if rhs==2 then
        sensors=[1:ny];
    end,
    if rhs==2 || rhs==3 then
        known=[];
    end,
    ts=sys(7);
    
    nsens = length(sensors); 
    nknown = length(known);
    [nl,ml] = size(l);
    if (ml~=nsens) then
        error(msprintf(gettext("%s:sensors length and l matrix no of columns should be same","estim")))
    end,
    if (nl~=nx) then
        error(msprintf(gettext("%s :A and l should have same number of rows"),"estim"))
    end
    b = b(:, known);
    c = c(sensors, :);
    d = d(sensors, known);
    inputs = [1:nsens] + nu; 
    states = [1:nx] + ny;
    //pause
    m = length (known);
    n = size(a,1);
    p = length (sensors);
    //state spaces of the estimator
    ae = a-l*c;
    be = [b-l*d,l];
    ce = [c;eye(nx,nx)];
    de = [d, zeros(p, p); zeros(n, m), zeros(n, p)];
    //dss(ae,be,ce,de,ts)
    
endfunction

