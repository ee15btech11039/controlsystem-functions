
function [N]=nmodels(sysarray)
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
        for i=1:nx
            for j=1:nu
                if and(typeof(sysarray{i,j})<>['rational','state-space']) then
                    error(msprintf(gettext(" %s:Wrong type of input model a state-space or a rational model expected \n"),"nmodels"))
                end,
            end,
        end,
        N=prod(size(sysarray));
        return;
    end
endfunction

