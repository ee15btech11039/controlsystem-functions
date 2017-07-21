
function [N]=nmodels(sysarray)
    //returns number of models in a model array
    //
    //Calling Seqence
    //N=nmodels(sysarray)
    //
    //Parameters
    //sysarray : lti system or a cell array of lti systems
    //N: no of models
    //
    //Description   
    //[rsys] = balred(sys,n) 
    // n = nmodels(sys)returns the number of models in an array of dynamic 
    //system models or static models.
    //
    //Author
    //Ayush Kumar
    
    [lhs,rhs]=argn(0)
    
    if rhs~=1 then
        error(msprintf(gettext(" %s:Wrong number of input argument,%d input argument expected.\n"),"nmodels",1))
    end,
    
    if ~(typeof(sysarray)=='ss' ||typeof(sysarray)=='state-space'||typeof(sysarray)=='ce')then
        error(msprintf(gettext(" %s:Wrong type of input argument,input argument should be an array of dynamic sytem models or static models.\n"),"nmodels",1))
    end
    
    if typeof(sysarray)=='ss' then
        k=sysarray.A
        sizearr=size(k)
        N=prod(sizearr(:,3:$))
        return;
    end
    
    if typeof(sysarray)=='state-space' then
        N=1;
        return;
    end
    
    if typeof(sysarray)=='ce' then
        [nx,nu]=size(sysarray);
        for i=1:prod(size(sysarray))
            if and(typeof(sysarray{i})<>['rational','state-space','ss']) then
                error(msprintf(gettext(" %s:Wrong type of input model a state-space or a rational model expected \n"),"nmodels"))
                
            end,
        end,
        count=0;
        for i=1:prod(size(sysarray))
            if typeof(sysarray{i})=="ss" then
                k=sysarray{i}.A
                sizearr=size(k)
                count=count+prod(sizearr(:,3:$))
            end
            if typeof(sysarray{i})=="state-space" || typeof(sysarray{i})=="rational" then
                count=count+1
            end
        end
        
        N=count;
        return;
    end
endfunction

