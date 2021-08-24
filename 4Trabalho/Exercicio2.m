function Exercicio2(file_name, c, k_gmres, maxit)
  
  load (file_name);
  A = Problem.A;
  n = rows(A);
  x_real = ones(n, 1);
  b = A*x_real;
  
  # Informações da matriz original
  k = -1;
  if (c)
    k = cond(A);
  endif
  nzeros = nnz(A);
  
  m_title = erase(erase(file_name, ".mat"), " ");
  print_name = strjoin({"Resultados", m_title, m_title}, "/");
  f_name = strjoin({print_name, "txt"}, ".");
  file = fopen(f_name, "w");
  fprintf(file, "Analise da matriz %s para uso de precondicionadores\n\n", m_title);
  fprintf(file, "Número de linhas da matriz: %d\n", n);
  fprintf(file, "Número de elementos não nulos: %d\n", nzeros);
    
  diagonal = "";
  if (diagonal_dominante(A))
    fprintf(file, "A matriz é diagonal dominante\n");
    diagonal = "s";
  else
    fprintf(file, "A matriz não é diagonal dominante\n");
    diagonal = "n";
  endif
  
  fprintf(file, "\nNúmero de condicionamento da matriz: %e\n\n", k);
  
  spy(A);
  print(strcat(print_name, "_spyA.png"));
  clf;
  
  # Informações da matriz resultante do precondicionamento ILU
  variacoes = ["ILU(0)"; "ILU"];
  reordenamento = ["com"; "sem"];
  
  figure(1)
  img_name =  strjoin({print_name, "png"}, ".");  
  xlabel ("iteração");
  ylabel ("log(er)");
  title(strjoin({"Gráfico do resíduo: ", m_title}))
  hold on

  # Sem precondicionadores
  
        tol = 10e-11;
        
        fprintf(file, "\n\nPara k= %d, %d número máximo de iterações e %e de tolerância\n", k_gmres, maxit, tol);
        tic;
        [x,flag,relres,iter,resvec] = gmres(A,b,k_gmres,tol,maxit);
        tempo = toc();
        norma = norm(x, inf);
        
        fprintf(file, "Tempo de espera: %d segundos\n", tempo);
        string_flag = "";
        if (flag == 0)
          string_flag = "convergencia atingida";
        elseif (flag == 1)
          string_flag = "número maximo de iterações atingido";
        elseif (flag == 3)
          string_flag = strjoin({"erro estagnado em", num2str(relres),"depois de", int2str( length(resvec) ), "iterações"});
        endif
        fprintf(file, "Flag de convergência: %d [%s]\n", flag, string_flag);
        fprintf(file, "Número de iterações: %d\n", length(resvec));
        fprintf(file, "Valor final do resíduo relativo: %e\n", relres);
        fprintf(file, "Norma do máximo da solução: %e\n", norma);
        
        fprintf(file, "Código latex da linha na tabela:\n\n");

        fprintf(file, "[-1em] \\bfseries %s & %d & %d & %d & %d & %e & %e & %d s \\\\ & & & & & & & & \\\\ [-1em] \\hline \\\\\n"
              , m_title, n, k_gmres, flag, length(resvec), relres, norma, tempo);
              
        figure(1)
        plot(log(resvec), "linewidth", 2, strjoin({";resíduo ", m_title, "sem precondicionamento;"}));
  
  # Com precondicionadores
    
  for i = 1:rows(variacoes)
    for j = 1:rows(reordenamento)
      fprintf(file, "\nPrecondicionamento %s %s reordenamento\n", erase(variacoes(i, :), " "), reordenamento(j, :));
      M = [];
      tic;
      try
        RP = [];
        if (strcmp("com", reordenamento(j, :)))
          perm = symrcm(A);
          I = speye(n,n);
          P = I(perm,:);
          RP = P*A*P';
          R_b = P*b;
        endif
        
        if (strcmp("ILU(0)", variacoes(i, :)))
          opts.type = "nofill";
          if (strcmp("com", reordenamento(j, :)))
            [M1, M2] = ilu(RP, opts);
          else
            [M1, M2] = ilu(A, opts);
          endif
          M = M1*M2;
        else
          opts.type = "crout";
          opts.droptol = 10e-4;
          if (strcmp("com", reordenamento(j, :)))
            [M1, M2] = ilu(RP, opts);
          else
            [M1, M2] = ilu(A, opts);
          endif
          M = M1*M2;
        endif
        tempo_precond = toc();
        
        nzeros_M = nnz(M);
        k_M = -1;
        if(c)
          k_M = cond(M);
        endif
        
        fprintf(file, "Número de elementos não nulos: %d\n", nzeros_M);
          
        diagonal = "";
        if (diagonal_dominante(M))
          fprintf(file, "A matriz é diagonal dominante\n");
          diagonal_M = "s";
        else
          fprintf(file, "A matriz não é diagonal dominante\n");
          diagonal_M = "n";
        endif
        
        fprintf(file, "\nNúmero de condicionamento da matriz: %e\n\n", k_M);
        
        figure(2)
        spy(M);
        print(erase(strcat(strjoin({print_name, "spyM", variacoes(i, :), reordenamento(j, :) }, "_"), ".png"), " "));
        clf;
        
        tol = 10e-11;
        
        fprintf(file, "\n\nPara k= %d, %d número máximo de iterações e %e de tolerância\n", k_gmres, maxit, tol);
        tic;
        if (strcmp("com", reordenamento(j, :)))
          [x,flag,relres,iter,resvec] = gmres(RP,R_b,k_gmres,tol,maxit,M);
        else
          [x,flag,relres,iter,resvec] = gmres(A,b,k_gmres,tol,maxit,M);
        endif
        tempo = toc();
        norma = norm(x, inf);
        
        fprintf(file, "Tempo de espera: %d segundos\n", tempo);
        string_flag = "";
        if (flag == 0)
          string_flag = "convergencia atingida";
        elseif (flag == 1)
          string_flag = "número maximo de iterações atingido";
        elseif (flag == 3)
          string_flag = strjoin({"erro estagnado em", num2str(relres),"depois de", int2str( length(resvec) ), "iterações"});
        endif
        fprintf(file, "Flag de convergência: %d [%s]\n", flag, string_flag);
        fprintf(file, "Número de iterações: %d\n", length(resvec));
        fprintf(file, "Valor final do resíduo relativo: %e\n", relres);
        fprintf(file, "Norma do máximo da solução: %e\n", norma);
        
        fprintf(file, "Código latex da linha na tabela:\n\n");

        fprintf(file, "[-1em] \\bfseries %s & %d & %d & %s & %e & %s & %s & %d & %s & %e & %d s \\\\ & & & & & & & & & & \\\\ [-1em] \\hline \\\\\n"
              , m_title, n, nzeros, diagonal, k, erase(variacoes(i, :), " "), reordenamento(j, :), nzeros_M, diagonal_M, k_M, tempo_precond);

        fprintf(file, "[-1em] \\bfseries %s & %d & %s & %s & %d & %s & %d & %d & %d & %d & %e & %e & %d s \\\\ & & & & & & & & & & & \\\\ [-1em] \\hline \\\\\n"
              , m_title, n, erase(variacoes(i, :), " "), reordenamento(j, :), nzeros_M, diagonal_M, k_M, k_gmres , flag, length(resvec), relres, norma, tempo);

              
        # Imprime imagem
        figure(1)
        plot(log(resvec), "linewidth", 2, strjoin({";resíduo ", m_title, erase(variacoes(i, :), " "), reordenamento(j, :) , "reordenamento;"}));
      catch
        msg = lasterror.message;
        fprintf(file, "\n\n%s\n\n", msg);
      end_try_catch
    endfor
  endfor
  
  fclose(file);
  
  print(img_name);
  hold off
  close;
  
  figure(2);
  close;
  
endfunction