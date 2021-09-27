diretorios = ["/Direto/";"/GMRES/"];

problema = "escoamento";

N = ["51x21"; "101x101";"201x201"];
tipos = [2,3];

for i = 1:rows(diretorios)
        for j = 1:rows(N)
            load(strjoin( { "Resultados", strrep(diretorios(i, :), " ", ""), problema, "_", strrep(N(j, :), " ", ""),   ".mat"}, "" ) );
            
            norma = norm(u, inf);
            string_flag = "";
            if (flag == 0)
              string_flag = "convergencia atingida";
            elseif (flag == 1)
              string_flag = "número maximo de iterações atingido";
            elseif (flag == 3)
              string_flag = strjoin({"erro estagnado em", num2str(relres),"depois de", int2str( length(resvec) ), "iterações"});
            endif
            
            f_name = strjoin( { "Resultados", strrep(diretorios(i, :), " ", ""), problema, "_", strrep(N(j, :), " ", ""),   ".txt"}, "" );
            file = fopen(f_name, "w");
            fprintf(file, "Número de linhas da matriz: %d\n", n*n);
            fprintf(file, "Discretização: %d\n", n);
            
            fprintf(file, "Tempo de espera: %d segundos\n", tempo);
            fprintf(file, "Flag de convergência: %d [%s]\n", flag, string_flag);
            fprintf(file, "Número de iterações: %d\n", length(resvec));
            fprintf(file, "Valor final do resíduo relativo: %e\n", relres);
            fprintf(file, "Norma do máximo da solução: %e\n", norma);
            fclose(file);
        endfor 
endfor