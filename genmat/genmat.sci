function [p] = genmat(a)
    //creating generalised matrices
    //
    //Calling Seqence
    //m=genmat(a)
    //
    //Parameters
    //a:realp object
    //m:generalized matrix
     
    //Description   
    //Generalized matrices (genmat) are matrices that depend on tunable parameters.
    //we can also use generalized matrices for building generalized LTI models that
    // represent control systems having a mixture of fixed and tunable components.
    //It can be used for parameter studies,how the system behaviour changes with 
    //change in parameter etc.
    // A generalized matrix is a structure having attributes blocks,SamplingGrid 
    //and doublem.where
    //Blocks : a structure which keep tracks of the realp objects which are used 
    //in creating the generalized matrix.
    //SamplingGrid : Sampling grid for model arrays
    //doublem : stores the current value of created generalized matrix 

    //   
    //Examples
    //a=realp('a',5);
    //b=realp('b',1);
    //m=a+b+1;
    //m.Blocks;m.doublem;
    //
    //Author
    //Ayush Kumar
    [lhs,rhs]=argn(0)
    if rhs<0 || rhs>1 then
        error(msprintf(gettext("%s:wrong number of input arguments %d argument expected"),"genmat",1))
    end,
    if rhs==0 then
        a=[];
    end
    if typeof(a)=="state-space"||typeof(a)=="rational"||typeof(a)=="polynomial" then
         error(msprintf(gettext("%s:lti model cannot be converted to genmat"),"genmat"))
    end,
    if typeof(a)=="realp" then
        p=struct("Blocks",a,"Sampling Grid",struct(),"doublem",a.value)
    end,
    if typeof(a)=="constant" then
        p=struct("Blocks",struct(),"Sampling Grid",struct(),"doublem",a)
    end
    
endfunction   
/////////////////////////////overloading functions///////////////////////////////
function [x]=%realp_a_s(a,b)  //addition using overloading
    x=struct("Blocks",struct(a.name,a),"SamplingGrid",struct(),"doublem",a.value+b)
endfunction
function x=%s_a_realp(a,b) // number+realp object
    x=struct("Blocks",struct(b.name,b),"SamplingGrid",struct(),"doublem",a+b.value)  
endfunction
function x=%realp_a_realp(a,b)  
     x=struct("Blocks",struct(a.name,a,b.name,b),"SamplingGrid",struct(),"doublem",a.value+b.value)
endfunction
function [x]=%st_a_realp(a,b)  //add here a is a structure and b realp
    x=a
    x.Blocks(b.name)=b
    x.doublem=a.doublem+b.value
endfunction
function [x]=%realp_a_st(a,b)  //add here b is a structure and a realp
    x=b
    x.Blocks(a.name)=a
    x.doublem=b.doublem+a.value
endfunction
function [x]=%st_a_s(a,b)  //add here a is a structure and b number
    x=a
    x.doublem=a.doublem+b
endfunction
function [x]=%s_a_st(a,b)  
    x=b
    x.doublem=a+b.doublem
endfunction



//function [x]=%realp_a_st(a,b)
  //  x=b
    //x.Blocks.(strinb.aane)=b  
 //x.doabuem=a.value+b.doublem
//endfunction


//--------------------------------------------------------------------------------
function [x]=%realp_s_s(a,b)  
    x=struct("Blocks",struct(a.name,a),"SamplingGrid",struct(),"doublem",a.value-b)
endfunction
function x=%s_s_realp(a,b)
    x=struct("Blocks",struct(b.name,b),"SamplingGrid",struct(),"doublem",a-b.value)  
endfunction
function x=%realp_s_realp(a,b)  
     x=struct("Blocks",struct(a.name,a,b.name,b),"SamplingGrid",struct(),"doublem",a.value-b.value)
endfunction
function [x]=%st_s_realp(a,b)  
    x=a
    x.Blocks(b.name)=b
    x.doublem=a.doublem-b.value
endfunction
function [x]=%realp_s_st(a,b) 
    x=b
    x.Blocks(a.name)=a
    x.doublem=b.doublem-a.value
endfunction
function [x]=%st_s_s(a,b)  
    x=a
    x.doublem=a.doublem-b
endfunction
function [x]=%s_s_st(a,b)  
    x=b
    x.doublem=a-b.doublem
endfunction

//---------------------------------------------------------------------------------
function [x]=%realp_m_s(a,b)  
    x=struct("Blocks",struct(a.name,a),"SamplingGrid",struct(),"doublem",a.value*b)
endfunction
function x=%s_m_realp(a,b)
    x=struct("Blocks",struct(b.name,b),"SamplingGrid",struct(),"doublem",a*b.value)  
endfunction
function x=%realp_m_realp(a,b)  
     x=struct("Blocks",struct(a.name,a,b.name,b),"SamplingGrid",struct(),"doublem",a.value*b.value)
endfunction
function [x]=%st_m_realp(a,b)  
    x=a
    x.Blocks(b.name)=b
    x.doublem=a.doublem*b.value
endfunction
function [x]=%realp_m_st(a,b) 
    x=b
    x.Blocks(a.name)=a
    x.doublem=b.doublem*a.value
endfunction
function [x]=%st_m_s(a,b)  
    x=a
    x.doublem=a.doublem*b
endfunction
function [x]=%s_m_st(a,b)  
    x=b
    x.doublem=a*b.doublem
endfunction


//----------------------------------------------------------------------------------

function [x]=%realp_r_s(a,b)  
    x=struct("Blocks",struct(a.name,a),"SamplingGrid",struct(),"doublem",a.value/b)
endfunction
function x=%s_r_realp(a,b)
    x=struct("Blocks",struct(b.name,b),"SamplingGrid",struct(),"doublem",a/b.value)  
endfunction
function x=%realp_r_realp(a,b)  
     x=struct("Blocks",struct(a.name,a,b.name,b),"SamplingGrid",struct(),"doublem",a.value/b.value)
endfunction
function [x]=%st_r_realp(a,b)  
    x=a
    x.Blocks(b.name)=b
    x.doublem=a.doublem/b.value
endfunction
function [x]=%realp_r_st(a,b) 
    x=b
    x.Blocks(a.name)=a
    x.doublem=b.doublem/a.value
endfunction
function [x]=%st_r_s(a,b)  
    x=a
    x.doublem=a.doublem/b
endfunction
function [x]=%s_r_st(a,b)  
    x=b
    x.doublem=a/b.doublem
endfunction

//-----------------------------------------------------------------------------------
function [x]=%realp_p_s(a,b)  
    x=struct("Blocks",struct(a.name,a),"SamplingGrid",struct(),"doublem",a.value^b)
endfunction
function x=%s_p_realp(a,b)
    x=struct("Blocks",struct(b.name,b),"SamplingGrid",struct(),"doublem",a^b.value)  
endfunction
function x=%realp_p_realp(a,b)  
     x=struct("Blocks",struct(a.name,a,b.name,b),"SamplingGrid",struct(),"doublem",a.value^b.value)
endfunction
function [x]=%st_p_realp(a,b)  
    x=a
    x.Blocks(b.name)=b
    x.doublem=a.doublem^b.value
endfunction
function [x]=%realp_p_st(a,b) 
    x=b
    x.Blocks(a.name)=a
    x.doublem=b.doublem^a.value
endfunction
function [x]=%st_p_s(a,b)  
    x=a
    x.doublem=a.doublem^b
endfunction
function [x]=%s_p_st(a,b)  
    x=b
    x.doublem=a^b.doublem
endfunction
///////////////////////////////////////////////////////////////////////////////////////
function [x]=%realp_l_s(a,b)  
    x=struct("Blocks",struct(a.name,a),"SamplingGrid",struct(),"doublem",a.value\b)
endfunction
function x=%s_l_realp(a,b)
    x=struct("Blocks",struct(b.name,b),"SamplingGrid",struct(),"doublem",a\b.value)  
endfunction
function x=%realp_l_realp(a,b)  
     x=struct("Blocks",struct(a.name,a,b.name,b),"SamplingGrid",struct(),"doublem",a.value\b.value)
endfunction
function [x]=%st_l_realp(a,b)  
    x=a
    x.Blocks(b.name)=b
    x.doublem=a.doublem\b.value
endfunction
function [x]=%realp_l_st(a,b) 
    x=b
    x.Blocks(a.name)=a
    x.doublem=b.doublem\a.value
endfunction
function [x]=%st_l_s(a,b)  
    x=a
    x.doublem=a.doublem\b
endfunction
function [x]=%s_l_st(a,b)  
    x=b
    x.doublem=a\b.doublem
endfunction

/////////////////////////////////////////////////////////////////////////////////////

//needs further deination of overloading functions like  %s_c_st,%s_f_st for the case
//of matrix operations like [a+1 0; b 1] where a,b realp objects.







 
       
