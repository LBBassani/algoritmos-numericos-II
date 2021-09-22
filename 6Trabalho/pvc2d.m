function [u] = pvc2d(a,b,c,d,n,m, app);

   % Discretização do domínio
   hx = (b-a)/(n-1);
   hy = (d-c)/(m-1);
   x  = linspace(a,b,n);
   y  = linspace(c,d,m);

   % Dados para caracterização do problema, funï¿½ï¿½es conhecidas e condições de contorno
   [bx,by,gamma,fun,kappa,left,gleft,right,gright,bottom,gbottom,top,gtop] = aplicacao (x,y,n,m, app);

   % Calculo dos coeficientes vindos do método das diferenças finitas
   [ai,bi,ci,di,ei]= coeficientes(hx,hy,kappa,bx,by,gamma,n,m);
   
   % Montagem inicial do sistema linear
   [A,f] = sistema_linear(ai,bi,ci,di,ei,fun,n,m);

   % Aplicação das condições de contorno
   [A,f] = condicoes_contorno (A,f,n,m,left,gleft,right,gright,bottom,gbottom,top,gtop,ai,bi,ci,di,ei,hx,hy,kappa);
 
   % Resolução do sistema
   u = A\f;

   % Plote do gráfico
   grafico_solucao (u,x,y,n,m)

endfunction
