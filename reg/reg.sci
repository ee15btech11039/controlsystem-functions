function [ac,bc,cc,dc]=reg(sys,k,l,sensors,known,controls)
    //(returns state spaces of a regulator given state-feedback and estimator gain.
    //
    //Calling Sequence
    //[ac,bc,cc,dc]=reg(sys,k,l)
    //[ac,bc,cc,dc]=reg(sys,k,l,sensors,controls)
    //
    //Parameters
    //sys:lti model (nx states,nu inputs,ny outputs)
    //k  :state feddback controller gain matrix (length(controls)-by-nu)
    //l  :estimator gain matrix (nx-by-length(sensors))
    //sensors : subset of output y which needs to be measured(default [1:ny],all 
    //the outputs measured)
    //known : known input indices (default=[])
    //controls: control input indices   (default=[1:nu])
    //[ae,be,ce,de]:state spaces of the estimator
    //note:length(known)+length(controls) shoult be less than or equal to nu.
    
    //Description   
    //[ac,bc,cc,dc] =reg(sys,K,L) returns the state spaces of  a dynamic 
    //regulator rsys given a state-space model  sys of the plant, a state-
    //feedback gain matrix K, and an estimator gain matrix L. The gains K  
    //and L are typically designed using pole  placement or  LQG techniques.
    //The function reg can be used for both continuous and discrete time 
    //plants.This syntax assumes that all inputs of sys are controls, and 
    //all outputs are measured. 
    //The regulator rsys is obtained by connecting the state-feedback law u = –Kx 
    //and the state estimator with gain matrix L which can be created by the state 
    //spaces ac,bc,cc,dc.This regulator should be connected to the plant using 
    //positive feedback.

    //For a plant:  dx/dt=Ax+Bu,                                                                             y=Cx+Du
    //The regulator is given by:                                                                _                     _
    //dx/dt = (A−LC−(B−LD)K)x + Ly
    //      _
    //u = -Kx
    //[ac,bc,cc,dc] = reg(sys,K,L,sensors,known,controls) handles more general 
    //regulation problems where: The plant inputs consist of controls u, known 
    //inputs ud, and stochastic inputs w.Only a subset y of the plant outputs is 
    //measured.
    //
    //
    
    //Author
    //Ayush Kumar
    [lhs,rhs]=argn(0);
    if typeof(sys)~="state-space" then
        error(msprintf(gettext("%s : wrong type of input arguments %d"),"reg",1))
    end,
    if rhs<3||rhs>6 then
        error(msprintf(gettext("%s : wrong number of input arguments %d"),"reg"))
    end,
    
    [a,b,c,d]=abcd(sys);
    [nx,na] = size(a);
    [ny,nu] = size(d);
    
    if rhs==3 then
        sensors = [1:ny];
        known = [];
        controls = [1:nu];
    end
    if rhs==4 then
        known = [];
        controls = [1:nu];
    end
    if rhs==5 then
        controls=[1:nu]
        controls(known)=[];
    end
    
    for i=1:length(sensors)
        if sensors(i)>ny then
            error(msprintf(gettext("%s : some index in argument %d is out of range"),"reg",4))
        end,
    end,
    for i=1:length(known)
        if controls(i)>nu then
            error(msprintf(gettext("%s : some index in argument %d is out of range"),"reg",5))
        end,
    end,
    for i=1:length(controls)
        if controls(i)>nu then
            error(msprintf(gettext("%s : some index in argument %d is out of range"),"reg",6))
        end,
    end,
    if length(known)+length(controls)>nu then
        error(msprintf(gettext("%s :length known and lengthcontrols vector should be less than or equal to inputs in the system"),"reg",nu))
    end
    nsens = length(sensors);
    nknown = length(known);
    nfb = length(controls);
   [nk,mk] = size(k);
   [nl,ml] = size(l);
   
    if (nk~=nfb) then
        error(msprintf(gettext("%s: no. of controls and size of k matrix should be same"),"reg")) 
    end
    if (mk~=nx) then
        error(msprintf(gettext("%s: no. of columns in k should be equal to no. of states in sys"),"reg")); 
    end,
    if (ml~=nsens) then
        error(msprintf(gettext("%s: no. of sensors should be same as size of l"),"reg")); 
    end,
    if (nl~=nx) then
         error(msprintf(gettext("%s: no. of rows in l should be same as no. of states in sys"),"reg"));
    end
   
    
    fdbk = [1:nfb] + ny; 
    inputs = [1:nsens] + nu;
    
    //----continuous LQG controller-----//
    ac = a;
    bc = [b,l];
    cc = [c;k];
    dc = [d, zeros(ny,nsens);zeros(nfb,nu), zeros(nfb,nsens)];
    
    //--sensors and internal control feedback loops--//
    
    S=[ac,bc;cc,dc]
    e=[fdbk,sensors];
    f=-[controls,inputs]
    outp=e;
    inp=abs(f),sgn=sign(f);
    nout=length(outp);nin=length(inp)
    [nx,na] = size(ac);
    [ny,nu]=size(dc);
    Bu = S(:,[nx+inp]); 
    Cy = S([nx+outp],:);
    if ~isempty(Cy)
        Cy(sgn==-1,:)=-Cy(sgn==-1,:);  // negative feedback
        E = eye(nout,nout) - Cy(:,[nx+inp]);    //E=(I-D22)
        Cy = E\Cy;
        Sc = S + Bu*Cy;   //Closed Loop
        
        //closed loop system
        ac = Sc(1:nx, 1:nx);
        bc = Sc(1:nx, nx+1:nx+nu);
        cc = Sc(nx+1:nx+ny, 1:nx);
        dc = Sc(nx+1:nx+ny, nx+1:nx+nu);
    else
        ac=a; bc=b; cc=c; dc=d;
    end
    //extracting subsystem//
    [nx,na]=size(ac);
    e=[known,inputs];f=fdbk;
    states=1:nx;
    ac=ac(states,states);
    bc=bc(states,e);
    cc=-cc(f,states);
    dc=dc(f,e)
    
endfunction
