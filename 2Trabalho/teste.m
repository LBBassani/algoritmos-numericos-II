function teste(file_name)
  
    for k = 1:rows(file_name)
      load(file_name(k, :));
      A = Problem.A;
      
      if diagonal_dominante(A)
        printf("%s é Diagonal Dominante\n", file_name(k, :));
      else
        printf("%s não é Diagonal Dominante\n", file_name(k, :));
      endif
      
    endfor
end 

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