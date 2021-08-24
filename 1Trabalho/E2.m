function E2(file_name)
  load (file_name);
  A = Problem.A;
  
  x_real = ones(rows(A), 1);
  b = A*x_real;
  
  img_name = erase(file_name, ".mat");
  m_title = erase( erase (img_name, "M100/" ), "M1000/");
  img_name = strcat(strcat(img_name, "/"), m_title);
  print_name = strcat(img_name, "_iterativo.txt");
  file = fopen(print_name, "w");
  fprintf(file, "Analise da matriz %s \n\n", m_title);
  
  if (diagonal_dominante(A))
    fprintf(file, "A matriz � diagonal dominante\n\n");
  else
    fprintf(file, "A matriz n�o � diagonal dominante\n\n");
  endif
  
  fclose(file);
  last_tol = 0.001;
  last_iter = 1000;
  
  file = fopen(print_name, "a");
  fprintf(file, "Jacobi tol = %e, max_iter = %d\n", last_tol, last_iter);
  [x_jacobi,er_jacobi,iter_jacobi] = jacobi(A, b, last_tol, last_iter, file);
  fclose(file);
  
  file = fopen(print_name, "a");
  fprintf(file, "Seidel tol = %e, max_iter = %d\n", last_tol, last_iter);
  [x_seidel,er_seidel,iter_seidel] = sor(A, b, last_tol, last_iter , 1, file);
  fclose(file);
  
  omega = [2:-0.2:0.2];
  col = columns(omega);
  best_sor = 0;
  best_omega = -1;
  best_er = 1000000;
  best_iter = 10000000;
  convergiu = 0;
  for i=1:col
    fprintf(stdout, "SOR w = %e, tol = %e, max_iter = %d\n", omega(i), last_tol, last_iter);
    [x_sor,er_sor,iter_sor] = sor(A,b,last_tol,last_iter,omega(i), stdout);
    if (not(convergiu))
      if (iter_sor < last_iter)
        convergiu = 1;
      elseif (er_sor(iter_sor) < best_er)
        best_er = er_sor(iter_sor);
        best_omega = omega(i);
        best_sor = er_sor;
      endif
     endif
     
     if(convergiu && best_iter > iter_sor)
        best_iter = iter_sor;
        best_er = er_sor(iter_sor);
        best_omega = omega(i);
        best_sor = er_sor;
     endif
    
  endfor

  best_omega = 1.8;
  file = fopen(print_name, "a");  
  fprintf(file, "SOR w = %e, tol %e, max_iter = %d\n", best_omega, last_tol, last_iter);
  [x_sor,er_sor,iter_sor] = sor(A,b,last_tol,last_iter,best_omega, file);
  
  fclose(file);
  
  figure("name", img_name);
  hold on;
  
  iter_jacobi = 1:iter_jacobi;
  plot(iter_jacobi, log(er_jacobi), ";Jacobi;");
  
  iter_seidel = 1:iter_seidel;
  plot(iter_seidel, log(er_seidel), ";Seidel;");
  
  iter_sor = 1:iter_sor;
  plot(iter_sor, log(er_sor), strcat(strcat("; SOR w = ", num2str(best_omega)), ";"));
  
  hold off;
  
  print(strcat(img_name, "_iterativo.png"));
  
  printf("Encerrado\n");
  
endfunction;

function [x,er,iter]=sor(A,b,tol,nmaxiter,w, file)
	   tic;
     [n,n]=size(A);
     iter = 1;
     er(1) = 1.0;  
     x0 = zeros(n,1);
     x = x0;
     while (er(iter) > tol )&&(iter < nmaxiter)

            for i=1:n
                soma = 0.0;
                for j = 1:(i-1)
                   soma = soma + A(i,j)*x(j);
                endfor
                for j = (i+1):n
		                soma = soma + A(i,j)*x0(j);
                endfor
	              x(i) = w*(b(i) - soma)/A(i,i) + (1-w)*x0(i);
            endfor
            iter = iter + 1;
	          er(iter) = norm(x-x0,inf)/norm(x,inf);            
            x0 = x;
     endwhile;
     normx = norm(x,inf);
     fprintf(file, "SOR - w=%f\n",w);
     fprintf(file, "Convergencia obtida apos %d iteracoes\n",iter);
     fprintf(file, "Norma do erro relativo = %e\n",er(iter));
	   fprintf(file, "Tempo de execu��o = %f s\n\n", toc());
endfunction;

function [t] = diagonal_dominante(A);
  [n,n]= size (A);
   for (i=1:n)
     soma =0.0;
     for (j=1:n)
       soma += abs(A(i,j));
     endfor
     soma = soma-abs(A(i,i));
     if (soma >= abs(A(i,i)))
       t=0;
       return;
     endif
   endfor
   t=1;

endfunction

function [x,er,iter]=jacobi(A,b,tol,nmaxiter, file)
		 tic;
     [n,n]=size(A);
     iter = 1;
     er(1) = 1.0;  
     x0 = zeros(n,1);
     x = x0;
     while ((er(iter) > tol )&&(iter < nmaxiter))

            for i=1:n
                soma = 0.0;
                for j = 1:(i-1)
                   soma = soma + A(i,j)*x0(j);
                endfor
                for j = (i+1):n
		                soma = soma + A(i,j)*x0(j);
                endfor
	              x(i) = (b(i) - soma)/A(i,i);
            endfor
            iter = iter + 1;
	          er(iter) = norm(x-x0,inf)/norm(x,inf);            
            x0 = x;	    
     endwhile;
     normx = norm(x,inf);
     fprintf(file, "Jacobi - Convergencia obtida apos %d iteracoes\n",iter);
     fprintf(file, "Norma do erro relativo = %e\n",er(iter));
		 fprintf(file, "Tempo de execu��o = %f s\n\n", toc());
endfunction;