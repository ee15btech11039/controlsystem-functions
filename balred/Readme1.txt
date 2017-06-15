function description :-

computes a reduced-order approximation of the input LTI system
syntax
[rsys]=balred(sys,n)
where sys can be in state-space or rational form,and n can be positive number or vector.

note:-in case n is a vector or matrix rsys retuned will be an array
of dimension (rows(n)*col(n),1)
whose elements can be viewed by the command rsys{k,1},where k varies from 
1 to no. of elements in n matrix,each rsys{k,1} will return the reduced order corresponding to the n matrix


algorithm:-here we are using balanced truncation algo to get approximated reduced system


ref: http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=261486&tag=1



