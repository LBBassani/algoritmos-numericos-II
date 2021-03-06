function Exercicio1(file_name, c)
  load (file_name);
  A = Problem.A;
  n = rows(A);
  x_real = ones(n, 1);
  b = A*x_real;
  
  # Informa??es da matriz original
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
  fprintf(file, "N?mero de linhas da matriz: %d\n", n);
  fprintf(file, "N?mero de elementos n?o nulos: %d\n", nzeros);
    
  diagonal = "";
  if (diagonal_dominante(A))
    fprintf(file, "A matriz ? diagonal dominante\n");
    diagonal = "s";
  else
    fprintf(file, "A matriz n?o ? diagonal dominante\n");
    diagonal = "n";
  endif
  
  fprintf(file, "\nN?mero de condicionamento da matriz: %e\n\n", k);
  
  if(c)
    spy(A);
    print(strcat(print_name, "_spyA.png"));
    close;
  endif
  # Informa??es da matriz resultante do precondicionamento ICC
  # Referencia: https://octave.org/doc/v5.2.0/Iterative-Techniques.html

  variacoes = ["ICC(0)"; "ICC"];
  reordenamento = ["com"; "sem"];
  
  figure(1)
  img_name =  strjoin({print_name, "png"}, ".");  
  xlabel ("itera??o");
  ylabel ("log(er)");
  title(strjoin({"Gr?fico do res?duo: ", m_title}))
  hold on
    
  for i = 1:rows(variacoes)
    for j = 1:rows(reordenamento)
      fprintf(file, "\nPrecondicionamento %s %s reordenamento\n", erase(variacoes(i, :), " "), reordenamento(j, :));
      M = [];
      tic;
      try
        if (strcmp("com", reordenamento(j, :)))
          perm = symrcm(A);
          I = speye(n,n);
          P = I(perm,:);
          RP = P*A*P';
          R_b = P*b;
        endif
        
        if (strcmp("ICC(0)", variacoes(i, :)))
          opts.type = "nofill";
          if (strcmp("com", reordenamento(j, :)))
            M1 = ichol(RP, opts);
          else
            M1 = ichol(A, opts);
          endif
          M2 = M1';
          M = M1*M2;
        else
          opts.type = "ict";
          opts.droptol = 10e-4;
          if (strcmp("com", reordenamento(j, :)))
            M1 = ichol(RP, opts);
          else
            M1 = ichol(A, opts);
          endif
          M2 = M1';
          M = M1*M2;
        endif
        tempo_precond = toc();
        
        nzeros_M = nnz(M);
        k_M = -1;
        if(c)
          k_M = cond(M);
        endif
        
        fprintf(file, "N?mero de elementos n?o nulos: %d\n", nzeros_M);
          
        diagonal = "";
        if (diagonal_dominante(M))
          fprintf(file, "A matriz ? diagonal dominante\n");
          diagonal_M = "s";
        else
          fprintf(file, "A matriz n?o ? diagonal dominante\n");
          diagonal_M = "n";
        endif
        
        fprintf(file, "\nN?mero de condicionamento da matriz: %e\n\n", k_M);
        
        if(c)
          figure(2)
          spy(M);
          print(erase(strcat(strjoin({print_name, "spyM", variacoes(i, :), reordenamento(j, :) }, "_"), ".png"), " "));
          close;
        endif      
        
        tol = 10e-11;
        maxit = 10e4;
        
        fprintf(file, "\n\nPara %d n?mero m?ximo de itera??es e %e de toler?ncia\n", maxit, tol);
        tic;
        if (strcmp("com", reordenamento(j, :)))
          [x,flag,relres,iter,resvec] = pcg(RP,R_b,tol,maxit, M);
        else
          [x,flag,relres,iter,resvec] = pcg(A,b,tol,maxit, M);
        endif
        tempo = toc();
        norma = norm(x, inf);
        
        fprintf(file, "Tempo de espera: %d segundos\n", tempo);
        string_flag = "";
        if (flag == 0)
          string_flag = "convergencia atingida";
        elseif (flag == 1)
          string_flag = "n?mero maximo de itera??es atingido";
        elseif (flag == 3)
          string_flag = strjoin({"erro estagnado em", num2str(relres),"depois de", int2str( length(resvec) ), "itera??es"});
        endif
        fprintf(file, "Flag de converg?ncia: %d [%s]\n", flag, string_flag);
        fprintf(file, "N?mero de itera??es: %d\n", length(resvec));
        fprintf(file, "Valor final do res?duo relativo: %e\n", relres);
        fprintf(file, "Norma do m?ximo da solu??o: %e\n", norma);
        
        fprintf(file, "C?digo latex da linha na tabela:\n\n");

        fprintf(file, "[-1em] \\bfseries %s & %d & %d & %s & %e & %s & %s & %d & %s & %e & %d s \\\\ & & & & & & & & & & \\\\ [-1em] \\hline \\\\\n"
              , m_title, n, nzeros, diagonal, k, erase(variacoes(i, :), " "), reordenamento(j, :), nzeros_M, diagonal_M, k_M, tempo_precond);

        fprintf(file, "[-1em] \\bfseries %s & %d & %s & %s & %d & %s & %d & %d & %e & %e & %e & %d s \\\\ & & & & & & & & & & & \\\\ [-1em] \\hline \\\\\n"
              , m_title, n, erase(variacoes(i, :), " "), reordenamento(j, :), nzeros_M, diagonal_M, k_M, flag, length(resvec), relres, norma, tempo);

              
        # Imprime imagem
        figure(1)
        plot(log(resvec), strjoin({";res?duo ", m_title, erase(variacoes(i, :), " "), reordenamento(j, :) , "reordenamento;"}));
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
endfunction