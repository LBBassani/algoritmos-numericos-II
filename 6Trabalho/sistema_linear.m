function [A,f] = sistema_linear (a,b,c,d,e,fun,n,m)
 
 N = n*m;
 A = sparse(N,N);
 f = zeros(N,1);


 for i = 1:N
   
   if i >= n + 1
      A(i,i-n) = d(i);
   endif
   if i > 1
      A(i,i-1) = b(i);
   endif
   
   A(i,i) = a(i);
   
   if i < N
      A(i,i+1) = c(i); 
   endif
   if i <= N - n
      A(i,i+n) = e(i);
   endif
   f(i) = fun(i); 
 endfor

 
endfunction

