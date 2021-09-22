function [u,x,y] = pvc2d(a,b,c,d,n,m, app);

   % Discretização do domínio
   hx = (b-a)/(n-1);
   hy = (d-c)/(m-1);
   x  = linspace(a,b,n);
   y  = linspace(c,d,m);

   % Dados para caracterização do problema, funï¿½ï¿½es conhecidas e condições de contorno
   [bx,by,gamma,fun,kappa,left,gleft,right,gright,bottom,gbottom,top,gtop] = app (x,y,n,m);

   % Calculo dos coeficientes vindos do método das diferenças finitas
   [ai,bi,ci,di,ei]= coeficientes(hx,hy,kappa,bx,by,gamma,n,m);
   
   % Montagem inicial do sistema linear
   [A,f] = sistema_linear(ai,bi,ci,di,ei,fun,n,m);

   % Aplicação das condições de contorno
   [A,f] = condicoes_contorno (A,f,n,m,left,gleft,right,gright,bottom,gbottom,top,gtop,ai,bi,ci,di,ei,hx,hy,kappa);
 
   % Resolução do sistema
   # u = A\f;
   k_gmres = floor(0.4*n*m); % estimar k cerca de 40% da dimensão da matriz
   tol = 10e-9;
   maxit = 10e4;
   [u,flag,relres,iter,resvec] = gmres(A,f,k_gmres,tol,maxit);

endfunction
