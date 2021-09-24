
function [A,f] = condicoes_contorno (A,f,n,m,left,gleft,right,gright,bottom,gbottom,top,gtop,a,b,c,d,e,hx,hy,kappa)
  
%
%  Contorno  "left" nos n�s I = 1,n+1,...,(m-1)*n+1 : 
% = 1 ent�o "gleft" � o valor prescrito "g" ;
% = 2 ent�o "gleft" � o valor do fluxo "sigma";
% = 3 ent�o "gleft" ser� um vetor ["alfa", "beta","gama"], onde gleft(1,1)=alfa; gleft(1,2)=beta; gleft(1,3)=gama;
%
%  Contorno  "right" nos n�s I = n,2*n,...,(m-1)*n, m*n:
% = 1 ent�o "gright" � o valor prescrito ;
% = 2 ent�o "gright" � o valor do fluxo sigma;
% = 3 ent�o "gright" ser� um vetor ["alfa", "beta","gama"], onde gright(1,1)=alfa; gright(1,2)=beta; gright(1,3)=gama;  
%
%  Contorno  "bottom" nos n�s I = 2,3,...,(n-1):
% = 1 ent�o "gbottom" � o valor prescrito ;
% = 2 ent�o "gbottom" � o valor do fluxo sigma;
% = 3 ent�o "gbottom" ser� um vetor ["alfa", "beta","gama"], onde gbottom(1,1)=alfa; gbottom(1,2)=beta; gbottom(1,3)=gama;  
%
%  Contorno  "top" nos n�s I = (m-1)*n+2,...,m*n-1:
% = 1 ent�o "gtop" � o valor prescrito ;
% = 2 ent�o "gtop" � o valor do fluxo sigma;
% = 3 ent�o "gtop" ser� um vetor ["alfa", "beta","gama"], onde gtop(1,1)=alfa; gtop(1,2)=beta; gtop(1,3)=gama;
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
       
       f(I) = f(I) + b(I).*(hx/kappa).*gleft(I);
       
       I(1) = [];
       a_index = sub2ind([N N], I, I-1);
       A(a_index) = b(I) = 0;
          
  case 3
       a_index = sub2ind([N N], I, I);
       alpha_index = sub2ind([N N], I, ones(1, length(I)));
       beta_index = sub2ind([N N], I, ones(1, length(I))*2);
       gamma_index = sub2ind([N N], I, ones(1, length(I))*3);
       
       A(a_index) = a(I) = a(I) + b(I).*( 1 - hx * (gleft(beta_index) ./ gleft(alpha_index) ) )';
       
       f(I) = f(I) - hx*( b(I).*(gleft( gamma_index) ./ gleft(alpha_index) )' );
       
       I(1) = [];
       a_index = sub2ind([N N], I, I-1);
       A(a_index) = b(I) = 0;
       
  otherwise
         printf("Erro na Condicao de contorno"); 
  end

 % Condicoes de contorno right
I = [n:n:m*n];
switch right  
case 1 
     [A, f] = valor_prescrito(A, f, gright, I);
case 2
     
     a_index = sub2ind([N N], I, I);
     A(a_index) = a(I) = a(I) + c(I);
       
     f(I) = f(I) + (hx/kappa)*(c(I).*gright(I));
     
     I(length(I)) = [];
     a_index = sub2ind([N N], I, I + 1);
     A(a_index) = c(I) = 0;
         
case 3
     a_index = sub2ind([N N], I, I);
     alpha_index = sub2ind([N N], I, ones(1, length(I)));
     beta_index = sub2ind([N N], I, ones(1, length(I))*2);
     gamma_index = sub2ind([N N], I, ones(1, length(I))*3);
     
     A(a_index) = a(I) = a(I) + c(I).*( 1 - hx * gright(beta_index) ./ gright(alpha_index) )';
     f(I) = f(I) - hx*( c(I).*(gright( gamma_index) ./ gright(alpha_index) )' );
       
     I(length(I)) = [];
     a_index = sub2ind([N N], I, I-1);
     A(a_index) = c(I) = 0;
otherwise
       printf("Erro na Condicao de contorno"); 
end

 % Condicoes de contorno bottom
I = [1:1:n];

% Garantindo valor prescrito nas quinas do dominio
if left == 1
  I(1) = [];
endif

if right == 1
  I(length(I)) = [];
endif

switch bottom
case 1
     [A, f] = valor_prescrito(A, f, gbottom, I); 
case 2
     
     a_index = sub2ind([N N], I, I);
     A(a_index) = a(I) = a(I) + d(I);
     
     f(I) = f(I) + (hy/kappa)*(d(I).*gbottom(I));
     
case 3
     a_index = sub2ind([N N], I, I);
     alpha_index = sub2ind([N N], I, ones(1, length(I)));
     beta_index = sub2ind([N N], I, ones(1, length(I))*2);
     gamma_index = sub2ind([N N], I, ones(1, length(I))*3);
     
     A(a_index) = a(I) = a(I) + d(I).*( 1 - hy * gbottom(beta_index) ./ gbottom(alpha_index) )';
     f(I) = f(I) - hy*( d(I).*(gbottom( gamma_index) ./ gbottom(alpha_index) )' );

otherwise
       printf("Erro na Condicao de contorno"); 
end

 % Condicoes de contorno top
I = [(m-1)*n+1:1:m*n];

% Garantindo valor prescrito nas quinas do dominio
if left == 1
  I(1) = [];
endif

if right == 1
  I(length(I)) = [];
endif

switch top
case 1 
     [A, f] = valor_prescrito(A, f, gtop, I);
case 2
     
     a_index = sub2ind([N N], I, I);
     A(a_index) = a(I) = a(I) + e(I);
     
     f(I) = f(I) + (hy/kappa)*(e(I).*gtop(I));
     
case 3
     
     a_index = sub2ind([N N], I, I);
     alpha_index = sub2ind([N N], I, ones(1, length(I)));
     beta_index = sub2ind([N N], I, ones(1, length(I))*2);
     gamma_index = sub2ind([N N], I, ones(1, length(I))*3);
     
     A(a_index) = a(I) = a(I) + e(I).*( 1 - hy * gtop(beta_index) ./ gtop(alpha_index) )';
     f(I) = f(I) - hy*( e(I).*(gtop( gamma_index) ./ gtop(alpha_index) )' );
     
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