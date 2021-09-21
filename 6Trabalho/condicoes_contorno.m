
function [A,f] = condicoes_contorno (A,f,n,m,left,gleft,right,gright,bottom,gbottom,top,gtop,a,b,c,d,e,hx,hy,kappa)
  
%
%  Contorno  "left" nos nós I = 1,n+1,...,(m-1)*n+1 : 
% = 1 então "gleft" é o valor prescrito "g" ;
% = 2 então "gleft" é o valor do fluxo "sigma";
% = 3 então "gleft" será um vetor ["alfa", "beta","gama"], onde gleft(1,1)=alfa; gleft(1,2)=beta; gleft(1,3)=gama;
%
%  Contorno  "right" nos nós I = n,2*n,...,(m-1)*n, m*n:
% = 1 então "gright" é o valor prescrito ;
% = 2 então "gright" é o valor do fluxo sigma;
% = 3 então "gright" será um vetor ["alfa", "beta","gama"], onde gright(1,1)=alfa; gright(1,2)=beta; gright(1,3)=gama;  
%
%  Contorno  "bottom" nos nós I = 2,3,...,(n-1):
% = 1 então "gbottom" é o valor prescrito ;
% = 2 então "gbottom" é o valor do fluxo sigma;
% = 3 então "gbottom" será um vetor ["alfa", "beta","gama"], onde gbotton(1,1)=alfa; gbotton(1,2)=beta; gbottom(1,3)=gama;  
%
%  Contorno  "top" nos nós I = (m-1)*n+2,...,m*n-1:
% = 1 então "gtop" é o valor prescrito ;
% = 2 então "gtop" é o valor do fluxo sigma;
% = 3 então "gtop" será um vetor ["alfa", "beta","gama"], onde gtop(1,1)=alfa; gtop(1,2)=beta; gtop(1,3)=gama;
%

N = n*m;

% Condicoes de contorno left
I = [1:n:(m-1)*n+1];

switch left
  case 1 
        [A, f] = valor_prescrito(A, f, gleft, I);
         
  case 2
       
       a_index = sub2ind([N N], I, I);
       A(a_index) = a(I) = a(I) + b(I);
       
       f(I) = f(I) + b(I)*(hx/kappa).*gleft(I);
       
       I(1) = [];
       a_index = sub2ind([N N], I, I-1);
       A(a_index) = b(I) = 0;
       
       A;
       f;
          
  case 3
         ...;
  otherwise
         printf("Erro na Condicao de contorno"); 
  end

 % Condicoes de contorno right
I = [n:n:m*n];
switch right  
case 1 
     [A, f] = valor_prescrito(A, f, gright, I);
case 2
     for i = I  
       if i < N
          A(i,i+1) = 0; 
       endif
         
         A(i,i) = a(i) + c(i);
         f(i) = f(i) + c(i)*(hx/kappa)*gright(i); 
     endfor    
case 3
      ...;
otherwise
       printf("Erro na Condicao de contorno"); 
end

 % Condicoes de contorno bottom
I = [1:1:n];
switch bottom
case 1
     [A, f] = valor_prescrito(A, f, gbottom, I); 
case 2
     for i = I  
       if i >= n + 1
          A(i,i-n) = 0;
       endif
         
         A(i,i) = a(i) + d(i);
         f(i) = f(i) + d(i)*(hy/kappa)*gbottom(i); 
     endfor
case 3
     ...;
otherwise
       printf("Erro na Condicao de contorno"); 
end

 % Condicoes de contorno top
I = [(m-1)*n+2:1:m*n-1];
switch top
case 1 
     [A, f] = valor_prescrito(A, f, gtop, I);
case 2
     for i = I  
       if i <= N - n
          A(i,i+n) = 0;
       endif
         
         A(i,i) = a(i) + e(i);
         f(i) = f(i) + e(i)*(hy/kappa)*gtop(i); 
     endfor
case 3
     ...;
otherwise
       printf("Erro na Condicao de contorno"); 
end

endfunction

function [A,f] = valor_prescrito (A, f, v, I)
    N = rows(A);
    a_index = sub2ind([N N], I, I);
    
    A(I, :) = 0; 
    A(a_index) = 1;
    f(I) = v;

endfunction