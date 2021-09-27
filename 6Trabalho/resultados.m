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
              string_flag = "n�mero maximo de itera��es atingido";
            elseif (flag == 3)
              string_flag = strjoin({"erro estagnado em", num2str(relres),"depois de", int2str( length(resvec) ), "itera��es"});
            endif
            
            f_name = strjoin( { "Resultados", strrep(diretorios(i, :), " ", ""), problema, "_", strrep(N(j, :), " ", ""),   ".txt"}, "" );
            file = fopen(f_name, "w");
            fprintf(file, "N�mero de linhas da matriz: %d\n", n*n);
            fprintf(file, "Discretiza��o: %d\n", n);
            
            fprintf(file, "Tempo de espera: %d segundos\n", tempo);
            fprintf(file, "Flag de converg�ncia: %d [%s]\n", flag, string_flag);
            fprintf(file, "N�mero de itera��es: %d\n", length(resvec));
            fprintf(file, "Valor final do res�duo relativo: %e\n", relres);
            fprintf(file, "Norma do m�ximo da solu��o: %e\n", norma);
            fclose(file);
        endfor 
endfor