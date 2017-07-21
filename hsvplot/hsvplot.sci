function [h]=hsvplot(sys)
    //Plot hankel singular values and return plot handles.
    //
    //Calling Seqence
    //[h]=hsvplot(sys)
    //hsvplot(sys)
    //
    //Parameters
    //sys   : a continuous or discrete time linear system.
    //h : plot handle

    //Description   
    //h=hsvplot(sys) plots the Hankel singular values of a given lti 
    //system and returns the plot handles h which can be used to further 
    //customize the plot.
 
    //Examples
    //sys=ssrand(1,1,8);t=0.1;
    //sys=dscr(sys,t);
    //h=hsvplot(sys);
 
    //Author
    //Ayush Kumar
    
    
    [lhs,rhs]=argn(0);
    //error checking ...............................................................
    if lhs<0 || lhs>1 then
        error(msprintf(gettext("%s: Wrong number of output arguments.\n"),"hsvplot"))
    end
    if rhs~=1  then
        error(msprintf(gettext("%s: Wrong number of input arguments:%d argument expected.\n"),"hsvplot",1))
    end,
    if (typeof(sys)=='rational') then
        if(degree(sys.num)>degree(sys.den))
            error(msprintf(gettext("The %s command cannot be used for models with more zeroes than poles","hankel")))
        end;
    end,

    if and(typeof(sys)<>['rational','state-space']) then
        error(msprintf(gettext("%s: Wrong type for input argument #%d: Linear state space or a transfer function expected.\n"),"hankel",1))
    end
    
    [ga,gs]=stabsep(sys);
    ninf=[];
    for i=1:size(ga.a,1)
        ninf=[ninf;%inf];
    end
    nk=hankel(sys); //squared hankel singular values
    nk=[ninf;nk]
    bar(nk^0.5)  //hankel singular values
    //note : in case system has unstable part infinite energies are shown 
    //as zeroes left to the first bar of stable modes.
    legend('stable modes');
    h=gca(); //plot handle
    
    h.x_label.text=_("State");
    h.y_label.text=_("State Energy");
    h.x_label.font_size=3;
    h.x_label.font_style=8;
    h.y_label.font_size=3;
    h.y_label.font_style=8;
    h.title.text=_("Hankel Singular Values");
    h.title.font_size=3;
    h.title.font_style=8;

    
endfunction

