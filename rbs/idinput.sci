//this function generates binary input signal depending upon the levels and 
//band
// N = [P nu M] gives a M*P-by-nu input, periodic with period P
//      and with M periods
//type:-'rbs'(random binary signal) here
//we are using a 8th order butterworth filter
//level:-binary outputs
//default P=1,NU=1,LEVELS=[-1 1],BAND=[0 0.5]
function[u]=idinput(N,types,band,levels)
    [lhs,rhs]=argn(0)
    if  rhs<4 then                    //checking the input given with function 
        levels=[-1 1]                 //arguments and allocating default value.if input is not  given
    end
    if  rhs<3 then
        band=[0 0.5]
    end
    if  rhs<2 then
        types='rbs'
    end
    if  size(N,2)==3 then
        P=N(1);
        nu=N(2);
        M=N(3);
    elseif  size(N,2)==2 then
        P=N(1);
        nu=N(2);
        M=1;
    elseif  size(N,2)==1 then
        P=N;
        nu=1;
        M=1;
    else
        error('erroininputargument');
    end

    if levels(2)<levels(1) then           //checking errors
        error('levels (1)should be less then levels(2)');
        end
    if type(types)~=10 then
        error('types can be rbs only input valid argument');
        end
    if(band(2)>0.5) then
        error('band should be in b/w 0-0.5');
    end
    
    function[vec]=randnu(NN,nuu)   //nuu channel random normally distributed signal generator
    [lhs,rhs]=argn(0)
    if  rhs<2 then                
        nuu=1                      
    end
    rand('seed',getdate('s'))
    for j=1:1:nuu
        for i= 1:NN
            vec(i,j)=rand(NN,'normal')
        end
    end
    endfunction


    if convstr(types)=='rbs' then
        u=randnu(5*P,nu);
        if ~and(band==[0 0.5]) then
            if(band(1)==0) then
                [hz]=iir(8,'lp','butt',[band(2) band(1)],[0 0]);  //8th order butterwoth filter
                num=hz(2);
                den=hz(3);
                for i=1:1:nu
                    y(:,i)=filter(num,den,u(:,i));
                end
            elseif(band(2)==0.5) then
                [hz]=iir(8,'hp','butt',[band(1) band(2)],[0 0]);  //8th order butterwoth filter
                num=hz(2);
                den=hz(3);
                for i=1:1:nu
                    y(:,i)=filter(num,den,u(:,i));
                end
            else
                [hz]=iir(8,'bp','butt',band,[0 0]);  //8th order butterwoth filter
                num=hz(2);
                den=hz(3);
                for i=1:1:nu
                    y(:,i)=filter(num,den,u(:,i));
                end
            end
            
            
            u = sign(y(2*P+1:3*P,:)); //taking the middle terms
            u = (levels(2)-levels(1))*(u+1)/2+levels(1);
        else
            u = sign(u(2*P+1:3*P,:));
            u = (levels(2)-levels(1))*(u+1)/2+levels(1);
            
        end
        
    else
        error('type can be(rbs)only')
    end
    if M>1 then                //generating periodic input if no. of periods>1
        uu = u;
        for i = 2:M
            u = [uu;u];
        end
    end
endfunction

    

    
    
