function E1(file_name)
  load (file_name);
  A = Problem.A;
  [L,U,P] = lu(A);
  
  img_name = erase(file_name, ".mat");
  m_title = erase( erase (img_name, "M100/" ), "M1000/");
  img_name = strcat(strcat(img_name, "/"), m_title);
  
  spy(A);
  print(strcat(img_name, "_spyA.png"));
  clf;
  spy(L);
  print(strcat(img_name, "_spyL.png"));
  clf;
  spy(U);
  print(strcat(img_name, "_spyU.png"));
  
  print_name = strcat(img_name, "_direto.txt");
  file = fopen(print_name, "w");
  fprintf(file, "Analise da matriz %s \n\n", m_title);
  
  preenchimento = 100 - 100*nnz(A)/( nnz(L) + nnz(U) );
  fprintf(file, "Taxa de preenchimento: %f%%\n", preenchimento );
  
  x_real = ones(rows(A), 1);
  b = A*x_real;
  x_barra = A\b;
  omega_x = x_real - x_barra;
  fprintf(file, "A distancia relativa entre a solução real e a aproximada é de %e\n", norm(omega_x, inf)/norm(x_barra, inf));
  
  omega_A = A - P*L*U;
  fprintf(file, "A distancia relativa entre a matriz original e a resultante da decomposição LU é de %e\n", norm( omega_A, inf)/norm( A, inf));
  
  omega_b = A*omega_x;
  fprintf(file, "A distancia relativa entre o vetor solução original e o resultante da decomposição LU é de %e\n", norm( omega_b, inf)/norm( b , inf));
  
  residuo = norm( b - A*x_barra, inf);
  fprintf(file, "A norma do resíduo é de %e\n", residuo);
  
  k = cond(A);
  fprintf(file, "O valor de condicionamento da matriz é %e\n", k);
  
  fclose(file);
  
endfunction;

