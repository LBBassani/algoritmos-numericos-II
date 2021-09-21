function [A,f] = sistema_linear (a,b,c,d,e,fun,n,m)
 
 N = n*m;
 A = sparse(N,N);
 f = zeros(N,1);


 for I = 1:N
   
   if I >= n + 1
      A(I,I-n) = d(I);
   endif
   if I > 1
      A(I,I-1) = b(I);
   endif
   
   A(I,I) = a(I);
   
   if I < N
      A(I,I+1) = c(I); 
   endif
   if I <= N - n
      A(I,I+n) = e(I);
   endif
   f(I) = fun(I); 
 endfor

 
endfunction

