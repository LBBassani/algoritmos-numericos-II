function [u] = pvc2d(a,b,c,d,n,m, app);

   % Discretiza��o do dom�nio
   hx = (b-a)/(n-1);
   hy = (d-c)/(m-1);
   x  = linspace(a,b,n);
   y  = linspace(c,d,m);

   % Dados para caracteriza��o do problema, fun��es conhecidas e condi��es de contorno
   [bx,by,gamma,fun,kappa,left,gleft,right,gright,bottom,gbottom,top,gtop] = aplicacao (x,y,n,m, app);

   % Calculo dos coeficientes vindos do m�todo das diferen�as finitas
   [ai,bi,ci,di,ei]= coeficientes(hx,hy,kappa,bx,by,gamma,n,m);
   
   % Montagem inicial do sistema linear
   [A,f] = sistema_linear(ai,bi,ci,di,ei,fun,n,m);

   % Aplica��o das condi��es de contorno
   [A,f] = condicoes_contorno (A,f,n,m,left,gleft,right,gright,bottom,gbottom,top,gtop,ai,bi,ci,di,ei,hx,hy,kappa);
 
   % Resolu��o do sistema
   u = A\f;

   % Plote do gr�fico
   grafico_solucao (u,x,y,n,m)

endfunction
