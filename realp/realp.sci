
function [var] = realp(paramname,initvalue,maximum,minimum,free)
    //real tunable parameter
    //
    //Calling Seqence
    //P=realp(name,initialvalue)
    //P=realp(name,initialvalue,max)
    //P=realp(name,initialvalue,max,min)
    //P=realp(name,initialvalue,max,min,free)
    //
    //Parameters
    //sys : siso or mimo lti system or system array(cell).
    //n : arraydimension
    //out:(cell)array of input systems
     
    //Description   
    //it creates a tunable real-valued parameter with name specified by 
    //name and initial value initvalue. it can be scalar- or matrix- valued.
    //max :  Upper bound for the parameter value.The dimension of the max 
    //property matches the dimension of the Value property.(default + inf)
    //min : Lower bound for the parameter value. The dimension of the min 
    //property matches the dimension of the Value property.(default - inf).
    //free:Boolean value specifying whether the parameter is free to be tuned.

    //   
    //Examples
    //a=realp('a',5);
    //
    //Author
    //Ayush Kumar
    [lhs,rhs]=argn(0)
    if rhs<2||rhs>5 then
        error(msprintf(gettext("%s : wrong number of input arguments \n","realp")));
    end,
    if typeof(paramname)~="string" then
        error(msprintf(gettext("%s: Wrong type of input argument :d a string expected \n"),"realp",1))
    end,
    if typeof(initvalue)~="constant" then
        error(msprintf(gettext("%s: Wrong type of input argument :d a constant expected \n"),"realp",2))
    end,
    if rhs==2 then
        maximum=ones(size(initvalue)(1),size(initvalue)(2))*%inf
        minimum=maximum*-1
        free=ones(size(initvalue)(1),size(initvalue)(2))
    end,
    if rhs==3 then
        minimum=-ones(size(initvalue)(1),size(initvalue)(2))*%inf
        free=ones(size(initvalue)(1),size(initvalue)(2))
    end,
    if rhs==4 then
        free=ones(size(initvalue)(1),size(initvalue)(2))
    end,
    
    if typeof(maximum)~="constant" then
        error(msprintf(gettext("%s: Wrong type of input argument :d a constant expected \n"),"realp",3))
    end,
    if typeof(minimum)~="constant" then
        error(msprintf(gettext("%s: Wrong type of input argument :d a constant expected \n"),"realp",3))
    end,
    if typeof(free)~="constant" then
        error(msprintf(gettext("%s: Wrong type of input argument :d a constant expected \n"),"realp",3))
    else
        if(free)~=0 then
            free=ones(size(initvalue)(1),size(initvalue)(2));
        end,
    end,
    
    var = tlist(["realp","name","value","max","min","free"],paramname,initvalue,maximum,minimum,free);//tlist declaration
endfunction

function %realp_p(var) //overloading to display the realp object
    f = fieldnames(var);
    if typeof(var(2))~="string" then
           error(msprintf(gettext("%s :name should  be of string type \n"),"realp"))
    end,
    if typeof(var(3))~="constant" then
           error(msprintf(gettext("%s :initial value should  be of constant type \n"),"realp"))
    end,
    if or(size(var(3))<>size(var(4)))||or(size(var(4))<>size(var(5)))||or(size(var(5))<>size(var(6)))||or(size(var(3))<>size(var(5)))||or(size(var(3))<>size(var(6)))||or(size(var(4))<>size(var(6))) then
        error(msprintf(gettext("%s: size of value,max,min,free cannot be changed"),"realp"))
    end,
    if typeof(var(4))~="constant" then
           error(msprintf(gettext("%s :max should  be of constant type \n"),"realp"))
    end,
    if typeof(var(5))~="constant" then
           error(msprintf(gettext("%s :min value should  be of constant type \n"),"realp"))
    end,
    if typeof(var(6))~="constant" then
           error(msprintf(gettext("%s :free should  be of constant type \n"),"realp"))
    end,
    [nx,nu]=size(var(3))
    if nx*nu>1 then
        disp("name :")
        disp(var(f(1)));
        disp("value:")
        disp(var(f(2)));
        disp("max :")
        disp(var(f(3)));
        disp("min :");
        disp(var(f(4)));
        disp("free :")
        disp(var(f(5)));
    else
        mprintf("     name :%s\n",var(2));
        mprintf("    value :%d\n",var(3));
        mprintf("      max :%d\n",var(4));
        mprintf("      min :%d\n",var(5));
        mprintf("     free :%d\n",var(6));
    end
    
endfunction    
    
    
    
    
    
    
    
    
    
    
    
