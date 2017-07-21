function [varargout]=rlocusplot(sys,kmax)
    //plots root locus and returns plot handle
    //
    //Calling Seqence
    //h=rlocusplot(H,Kmax)
    //
    //Parameters
    //H :- Siso linear system given by a transfer function or a state space 
    //representation.
    //Kmax real(maximum gain desired for the plot)
    //H: plot handle,which can be used to further customize the plot.
    //
    //Description   
    // Gives the root locus for a SISO linear system in state-space or transfer
    // form H(s) . This is the locus of the roots of 1+k*H(s)=1+k*N(s)/D(s), in 
    //the complex plane. For a selected sample of gains k <= kmax ,the imaginary
    // part of the roots of D(s)+k*N(s) is plotted vs the real part. 
    //   
    //Examples
    //sys=ssrand(1,1,5);
    //scf();
    //h=rlocusplot(sys);
    //h.font_size=2;
    //
    //Author
    //Ayush Kumar
    
    
    [lhs,rhs]=argn(0);
    if rhs <= 0 then   // demonstration
        n=real(poly([0.1-%i 0.1+%i,-10],"s"));
        d=real(poly([-1 -2 -%i %i],"s"));
        kmax=80;
        evans(n,d,kmax);
        if lhs==1 then
            varargout(1)=gca();
        end,
        return
    end,
    
    if rhs>2 then
        error(msprintf(gettext("%s: wrong number of input arguments atmost %d  arguments expected","rlocusplot",2)))
    end,
    if lhs>1 then
        error(msprintf(gettext("%s: wrong number of output arguments","rlocusplot")))
    end

    //numerator and denominator
    select typeof(sys)
    case "polynomial"  then
        n=sys
        d=1
        if rhs==1 then 
            kmax=0,end
    case "rational" then
        [n,d]=sys(2:3)
        if rhs==1 then 
            kmax=0,
        end
    case "state-space" then
        sys=ss2tf(sys);
        [n,d]=sys(2:3)
        if rhs==1 then kmax=0,end
        n=clean(n);d=clean(d);
    else
        error(msprintf(gettext("%s: wrong type of input argument %d a real dynamic system or polynomial expected"),"rlocusplot",1));
    end
    if prod(size(n))<>1 then
        error(msprintf(_("%s: Wrong value for input argument #%d: Single input, single output system expected.\n"),"rlocusplot",1));
    end
    if degree(n)==0&degree(d)==0 then
        error(msprintf(_("%s: The given system has no poles and no zeros.\n"),"rlocusplot"));
    end
    //plotting the rlocusplot and returning plot handle
    
    evans(n,d,kmax)
        
    h=gca();
    h.title.text=_("Root locus plot");
    h.title.font_style=8;
    h.title.font_size=3;
    h.x_label.font_size=3;
    h.y_label.font_size=3;
    h.x_label.font_style=8;
    h.y_label.font_style=8;
    if lhs==1 then
        varargout(1)=h;
    end,
endfunction

