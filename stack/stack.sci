function [out]=stack(n,varargin)
    //build model (cell)array by stacking models along array dimensions
    //
    //Calling Seqence
    //out=stack(n,sys1,sys2.....)
    //
    //Parameters
    //sys : siso or mimo lti system or system array(cell).
    //n : arraydimension
    //out:(cell)array of input systems
     
    //Description   
    //out = stack(n,sys1,sys2,...) produces an array of dynamic system models
    // out by stacking (concatenating) the models  (or arrays) sys1,sys2,... 
    //along the array dimension. All system models must have the same number of inputs
    // and outputs (the same I/O dimensions)
    //   
    //Examples
    //s=%s;
    //sys1=s/(s*s+5s+1)
    //sys2=s/(s+1)
    //out1=stack(1,sys1,sys2)
    //out1{:,:}
    //out2=stack(3,sys1,sys2)
    //out2{:,:}
    //
    //Author
    //Ayush Kumar
    
    ///////////////////////////////////////////////////////////////////////////////////////
    [lhs,rhs]=argn(0)
    for i=1:rhs-1
        if (typeof(varargin(i))~="state-space"&&typeof(varargin(i))~="rational"&&typeof(varargin(i))~="polynomial"&&typeof(varargin(i))~="ce") then 
            error(msprintf(gettext("%s :input argument should be a lti sysytem or lti system array"),"stack"))
        end,
    end,
    
    flag3=0
    for i=1:rhs-1
        if typeof(varargin(i))=="ce" then flag3=1 end,
    end
    
    if flag3==1 then
        for i=1:rhs-1 //in case of a cell array input ,all inputs should be a cell array
            if typeof(varargin(i))~="ce" then
                error(msprintf(gettext("%s:wrong type of input argument input"),"stack"))
            end
            
        end
        
    end,
    
    
    count=0;
    for i=1:rhs-1
        if typeof(varargin(i))=="rational" || typeof(varargin(i))=="polynomial" then
            count=count+1;
        end,
        if typeof(varargin(i))=="state-space" then
            count=count+2;
        end,
        if typeof(varargin(i))=="ce" then
            count=count+3;
        end,
    end,
    /////pure rational or polynomial input models
    if count==rhs-1 then
        for i=1:rhs-1
            testmat{i,1}=varargin(i); //storing data in a cell
        end
    end
    
    /////////mixture state-space and rational input models///////////////
    
    
    if count>(rhs-1) then
        for i=1:rhs-1
            //if flag==1 then
                if typeof(varargin(i))=="rational" || typeof(varargin(i))=="polynomial" then
                    testmat{i,1}=tf2ss(varargin(i)) //converting to state-space
                elseif typeof(varargin(i))=="state-space" then
                    testmat{i,1}=varargin(i);
                elseif typeof(varargin(i))=="ce" then
                    tempcell=varargin(i);
                    for j=1:prod(size(tempcell))
                        if typeof(tempcell{j})=="rational" || typeof(tempcell{j})=="polynomial" then
                            k=tf2ss(tempcell{j})
                            tempcell{j}=k
                        end
                    end
                    testmat{i,1}=tempcell;
                   
                end
                
        end 
        //k=testmat{1,1}
        //[nxip,nuip]=size(k{1,1})
        //[nxop,nuop]=size(testmat{1,1}.C)
        for i=1:rhs-1
            if typeof(testmat{i,1})=="state-space" then 
                [nxip,nuip]=size(testmat{1,1})
                if or(size(testmat{i,1})<>[nxip,nuip])   then       
                    error(msprintf(gettext("%s:input systems should have same input output matrix sizes"),"stack"))
                end
            end
            if typeof(testmat{i,1})=="ce" then 
                //k=testmat{1,1}
                //[nxip,nuip]=size(k{1,1})
                for j=1:prod(size(testmat{i,1}))-1
                    tempcell=testmat{i,1};
                    if or(size(tempcell{j})<>size(tempcell{j+1}))   then       
                        error(msprintf(gettext("%s:input systems should have same input output matrix sizes"),"stack"))
                    end,
                end,
            end,
        end,
    end,
    
    mat=[]
    for i=1:n-1
        mat(i)=1
    end,
    mat(n)=rhs-1;
    if n==1 then
        tempmat=cell(mat,n);
    else
        tempmat=cell(mat);
    end
    if n==1 then
        for i=1:rhs-1
            for j=1:n   //adding system in tempmat along array dimension 
                tempmat{i,j}=testmat{i+j-1,1}
            end,
        end,
    else
        for i=1:rhs-1
            tempmat{$,i}=testmat{i,1};
        end,
    end
    
    out=tempmat;
    
endfunction
//out{:,:}
