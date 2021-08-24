files = {"oscil_dcop_02"};

for i = 1:1
    matriz = files{i};
    titulo = ["matriz " matriz]

    load([matriz ".mat"])

    A = Problem.A;

    n = rows(A)
    ordem = floor(log10(n));
    magnitude = ["10^" num2str(ordem)];
    naonulos = nnz(A)

    b = A * ones(n, 1);
    tolordem_array = [6 9 12];
    k_array = [10 100 430];
    maxit = 1000

    figure("name", matriz)
    set (0, "defaultaxesfontname", "Helvetica")
    xlabel ("Iteracoes");
    ylabel ("log(resvec)");
    hold on
    for ki = 1:3
        k = k_array(ki)

        for toli = 1:3
            rtol = 1 * (10^(-(tolordem_array(toli))))
            tic
            [x, flag, relres, iter, resvec] = gmres(A, b, k, rtol, maxit);
            toc
            flag

            [iter2, _] = size(resvec)
            printf("Residuo minimo Generalizado - Convergencia obtida apos %d ciclos e %d iteracoes\n", iter(1, 1), iter2);

            printf("Norma da solucao = %f\n", norm(x, inf));
            
            plot([1:iter2], log(resvec), strcat(";residuo:", num2str(k), "-", num2str(tolordem_array(toli)), ";"));

            printf("\n\n")
        endfor

    endfor
    t = title(["Historico de Convergencia da " titulo]);
    set(t, 'Interpreter', 'none')
    set(0, 'DefaultTextInterpreter', 'none')
    saveas(1, strcat(matriz, ".jpg"));
    hold off
endfor
