function Exercicio()
   I = [1,2,3];
   N = [100];
   for i = I
     for n_m = N
       
         % Análise assintotica com n = m
         a = c = 0;
         b = d = 1;
         [u,flag,relres,iter,resvec,x,y] = pvc2d(a,b,c,d,n_m,n_m, @(x,y,n,m) analise_assintotica(x,y,n,m,i));
         
         figure()
         grafico_solucao(u,x,y,n_m,n_m)
         hgsave (strjoin({"Figuras/Direto/assintotica_"; int2str(n_m); "_tipo"; int2str(i); ".ofig"}, ""))
         close()
         save(strjoin({"Resultados/Direto/assintotica_"; int2str(n_m); "_tipo"; int2str(i); ".mat"}, ""), "u","flag","relres","iter","resvec");
         
      endfor
     endfor
     
     L = W = 1;
     
     for i = I
       for n_m = N
         
         a = c = 0;
         b = L;
         d = W;
         
         [u,flag,relres,iter,resvec,x,y] = pvc2d(a,b,c,d,n_m,n_m, @(x,y,n,m) resfriador_bidimensional(x,y,n,m,i));
       
         figure()
         grafico_solucao(u,x,y,n_m,n_m)
         hgsave (strjoin({"Figuras/Direto/resfriador_"; int2str(n_m); "_tipo"; int2str(i); ".ofig"}, ""))
         close()
         save(strjoin({"Resultados/Direto/resfriador_"; int2str(n_m); "_tipo"; int2str(i); ".mat"}, ""), "u","flag","relres","iter","resvec");
         
       endfor
     endfor
     
endfunction

function [bx,by,gamma,fun,kappa,left,gleft,right,gright,bottom,gbottom,top,gtop] = resfriador_bidimensional(x,y,n,m,ver)
    
    N = n*m;
    T = 2;
    kappa = 1.0;
    c = 1;
    u_ref = 70;
    
    I = [1:1:N];
    bx = sparse(N, 1);
    by = sparse(N, 1);
    gamma = sparse(N, 1);
    fun = sparse(N, 1);
    
    gamma(I) = 2*c/T;
    fun(I) = 2*c*u_ref/T;
    
    switch ver
       case 1
           
         left = 1;
         gleft = 200.0;

         right = 1;
         gright = 70.0;

         bottom = 1; 
         gbottom = 70.0;

         top = 1;
         gtop = 70.0;
       case 2
         
         left = 1;
         gleft = 200;
         
         right = 1;
         gright = 70;

         bottom = 2; 
         gbottom = sparse(n*m, 1);

         top = 2;
         gtop = sparse(n*m, 1);
         
       case 3
           
         left = 3;
         I = [1:n:(m-1)*n+1];
         gleft = sparse(n*m, 3);
         
         alpha_index = sub2ind([N N], I, ones(1, length(I)));
         beta_index = sub2ind([N N], I, ones(1, length(I))*2);
         gamma_index = sub2ind([N N], I, ones(1, length(I))*3);
         
         gleft(alpha_index) = kappa;
         gleft(beta_index) = c;
         gleft(gamma_index) = c*u_ref;

         right = 1;
         gright = 70.0;

         bottom = 1; 
         gbottom = 70.0;

         top = 1;
         gtop = 70.0;
     endswitch
    
endfunction

function [bx,by,gamma,fun,kappa,left,gleft,right,gright,bottom,gbottom,top,gtop] = analise_assintotica(x,y,n,m,ver)
    
    N = n*m;
    [X, Y] = meshgrid(x, y);
    x_index = @(I) rem(I - 1, n) + 1;
    y_index = @(I) floor((I - 1)/n) + 1;
    
    u = @(x, y) 10 .* (x .* y) .* (1 - x) .* (1 - y) .* e.^(x.^4.5);
    u_x = @(x, y) -45 .* e.^(x.^4.5) .* (0.222222 - 0.444444 .* x + x.^4.5 - x.^5.5) .* (-1 + y) .* y;
    u_x_x = @(x, y) 337.5 .* e.^(x.^4.5) .* (0.0592592 - 0.733333 .* x.^3.5 + x.^4.5 - 0.6 .* x.^8 + 0.6 .* x.^9) .* (-1 + y) .* y;
    u_y = @(x, y) 10 .* e.^(x.^4.5) .* (-1 + x) .* x .* (-1 + 2 .* y);
    u_y_y = @(x, y) 20 .* e.^(x.^4.5) .* (-1 + x) .* x;
    
    kappa = 1.0;
    I = [1:1:N];
    bx = sparse(N, 1);
    by = sparse(N, 1);
    gamma = sparse(N, 1);
    
    bx(I) = 1;
    by(I) = 20*y(y_index(I));
    gamma(I) = 1;
    
    ureal = u(X, Y);
    fun = -kappa * (u_x_x(X, Y)'(:) + u_y_y(X, Y)'(:)) + bx .* u_x(X, Y)'(:) + by .* u_y(X, Y)'(:) + ureal'(:);

    switch ver
       case 1 # valor prescrito
           
         left = 1;
         gleft = 0.0;

         right = 1;
         gright = 0.0;

         bottom = 1; 
         gbottom = 0.0;

         top = 1;
         gtop = 0.0;
       
        case 2 # fluxo prescrito
         
         # fluxo prescrito a esquerda
         left = 2;
         I = [1:n:(m-1)*n+1];
         gleft = sparse(n*m,1);
         
         x_i = x(x_index(I));
         y_i = y(y_index(I));
         
         du = u_x(x_i, y_i);
         gleft(I) = du*kappa;
         
         # fluxo prescrito a direita
         right = 2;     
         I = [n:n:m*n];
         gright = sparse(n*m,1);

         x_i = x(x_index(I));
         y_i = y(y_index(I));
         
         du = u_x(x_i, y_i);     
         gright(I) = -du*kappa;
         
         # fluxo prescrito para baixo
         bottom = 2;
         I = [2:1:n-1];
         gbottom = sparse(n*m, 1);
         
         x_i = x(x_index(I));
         y_i = y(y_index(I));
         
         du = u_y(x_i, y_i);
         gbottom(I) = du./kappa;

         # fluxo prescrito para cima
         top = 2;
         I = [(m-1)*n+2:1:m*n-1];
         gtop = sparse(n*m, 1);
         
         x_i = x(x_index(I));
         y_i = y(y_index(I));
         
         du = u_y(x_i, y_i);
         gbottom(I) = -du./kappa;
         
        case 3 # misto
          
         # esquerda 
         left = 3;
         I = [1:n:(m-1)*n+1];
         gleft = sparse(n*m, 3);
         
         alpha_index = sub2ind([N N], I, ones(1, length(I)));
         beta_index = sub2ind([N N], I, ones(1, length(I))*2);
         gamma_index = sub2ind([N N], I, ones(1, length(I))*3);
         
         gleft(alpha_index) = 1.0;
         gleft(beta_index) = 1.0;
         
         x_i = x(x_index(I));
         y_i = y(y_index(I));
         gleft(gamma_index) = -u_x(x_i, y_i) + u(x_i, y_i);

         # direita
         right = 3;
         I = [n:n:m*n];
         gright = sparse(n*m, 3);

         x_i = x(x_index(I));
         y_i = y(y_index(I));
         
         alpha_index = sub2ind([N N], I, ones(1, length(I)));
         beta_index = sub2ind([N N], I, ones(1, length(I))*2);
         gamma_index = sub2ind([N N], I, ones(1, length(I))*3);
         
         gright(alpha_index) = 1.0;
         gright(beta_index) = 1.0;
         gright(gamma_index) = u_x(x_i, y_i) + u(x_i, y_i);

         # para baixo
         bottom = 3; 
         I = [2:1:n-1];
         gbottom = sparse(n*m, 3);
         
         x_i = x(x_index(I));
         y_i = y(y_index(I));
         
         alpha_index = sub2ind([N N], I, ones(1, length(I)));
         beta_index = sub2ind([N N], I, ones(1, length(I))*2);
         gamma_index = sub2ind([N N], I, ones(1, length(I))*3);
         
         gbottom(alpha_index) = 1.0;
         gbottom(beta_index) = 1.0;
         gbottom(gamma_index) = -u_y(x_i, y_i) + u(x_i, y_i);

         # para cima
         top = 3;
         I = [(m-1)*n+2:1:m*n-1];
         gtop = sparse(n*m, 3);
         
         x_i = x(x_index(I));
         y_i = y(y_index(I));
         
         alpha_index = sub2ind([N N], I, ones(1, length(I)));
         beta_index = sub2ind([N N], I, ones(1, length(I))*2);
         gamma_index = sub2ind([N N], I, ones(1, length(I))*3);
         
         gtop(alpha_index) = 1.0;
         gtop(beta_index) = 1.0;
         gtop(gamma_index) = u_y(x_i, y_i) + u(x_i, y_i);
   
 endswitch 
endfunction