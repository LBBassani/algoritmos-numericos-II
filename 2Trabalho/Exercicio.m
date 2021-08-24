function Exercicio(file_name, c)
  load (file_name);
  A = Problem.A;
  n = rows(A);
  x_real = ones(n, 1);
  b = A*x_real;
  
  
  k = -1;
  if (c)
    k = cond(A);
  endif
  nzeros = nnz(A);
  
  m_title = erase(file_name, ".mat");
  print_name = strjoin({"Resultados", m_title, m_title}, "/");
  f_name = strjoin({print_name, "txt"}, ".");
  file = fopen(f_name, "w");
  fprintf(file, "Analise da matriz %s na resolução pelo método dos Gradientes Conjugados\n\n", m_title);
  fprintf(file, "Número de linhas da matriz: %d\n", n);
  fprintf(file, "Número de elementos não nulos: %d\n", nzeros);
  fprintf(file, "\nNúmero de condicionamento da matriz: %e\n", k);
  
  tol = [10e-6; 10e-11];
  maxit = [10000; 5000];
  
  for i = 1:rows(maxit)
    for j = 1:rows(tol)
      fprintf(file, "\n\nPara %d número máximo de iterações e %e de tolerância\n", maxit(i), tol(j));
      tic;
      [x,flag,relres,iter,resvec] = pcg(A,b,tol(j),maxit(i));
      tempo = toc();
      norma = norm(x, inf);
      
      fprintf(file, "Tempo de espera: %d segundos\n", tempo);
      string = "";
      if (flag == 0)
        string = "convergencia atingida";
      elseif (flag == 1)
        string = "número maximo de iterações atingido";
      elseif (flag == 3)
        string = strjoin({"erro estagnado em", num2str(relres),"depois de", int2str( length(resvec) ), "iterações"});
      endif
      fprintf(file, "Flag de convergência: %d [%s]\n", flag, string);
      fprintf(file, "Número de iterações: %d\n", length(resvec));
      fprintf(file, "Valor final do resíduo relativo: %e\n", relres);
      fprintf(file, "Norma do máximo da solução: %e\n", norma);
      
      fprintf(file, "Código latex da linha na tabela:\n\n");
      fprintf(file, "[-1em] \\bfseries %s & %d & %d & %d & %d & %e & %e & %e & %d s \\\\ & & & & & & & & \\\\ [-1em] \\hline \\\\\n"
            , m_title, n, nzeros, flag, length(resvec), relres, norma, k, tempo);
      
      
      img_name = strjoin({print_name, int2str(maxit(i)), int2str(log10(tol(j)) - 1)}, "_");
      img_name = strcat(img_name, ".png");
      plot(log(resvec), ";resíduo;");
      xlabel ("iteração");
      ylabel ("log(er)");
      title (strjoin({"Gráfico do resíduo: ", m_title, " \nmaxit =", int2str(maxit(i)), "tol =", num2str(tol(j))}));
      print(img_name);
      close;
    endfor
  endfor
  
  fclose(file);
  
end