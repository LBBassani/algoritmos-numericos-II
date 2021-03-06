function [u,flag,relres,iter,resvec,x,y] = pvc2d(a,b,c,d,n,m, app, direto);

   % Discretiza??o do dom?nio
   hx = (b-a)/(n-1);
   hy = (d-c)/(m-1);
   x  = linspace(a,b,n);
   y  = linspace(c,d,m);

   % Dados para caracteriza??o do problema, fun��es conhecidas e condi??es de contorno
   [bx,by,gamma,fun,kappa,left,gleft,right,gright,bottom,gbottom,top,gtop] = app (x,y,n,m);

   % Calculo dos coeficientes vindos do m?todo das diferen?as finitas
   [ai,bi,ci,di,ei]= coeficientes(hx,hy,kappa,bx,by,gamma,n,m);
   
   % Montagem inicial do sistema linear
   [A,f] = sistema_linear(ai,bi,ci,di,ei,fun,n,m);

   % Aplica??o das condi??es de contorno
   [A,f] = condicoes_contorno (A,f,n,m,left,gleft,right,gright,bottom,gbottom,top,gtop,ai,bi,ci,di,ei,hx,hy,kappa);
 
   % Resolu??o do sistema
   if direto
     u = A\f;
     flag = relres = iter = resvec = 0;
   else
     k_gmres = round(n*m*0.005); % estimar k cerca de 0.5% da dimens?o da matriz
     tol = 10e-9;
     maxit = 10e4;
     [u,flag,relres,iter,resvec] = gmres(A,f,k_gmres,tol,maxit);
   endif
   
endfunction
