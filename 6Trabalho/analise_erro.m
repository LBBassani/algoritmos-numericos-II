diretorios = ["/Direto/";"/GMRES/"];
problema = "assintotica";
N = [32,16,8,4];
tipos = [1,2,3];

a = c = 0;
b = d = 1;

for i = 1:rows(diretorios)
    for t = tipos
        for j = 1:length(N)
            load(strjoin( { "Resultados", strrep(diretorios(i, :), " ", ""), problema, "_", int2str(N(j)), "_tipo", int2str(t), ".mat"}, "" ) );
            
            % Valores reais
            f = @(x, y) 10 .* (x .* y) .* (1 - x) .* (1 - y) .* e.^(x.^4.5);
            x_index = @(I) rem(I - 1, n) + 1;
            y_index = @(I) floor((I - 1)/n) + 1;
            I = [1:N(j)];
            ureal = f(x(x_index(I)), y(y_index(I)));
            err = ureal - u;
            norma_err(j) = norm(err, Inf);
            h(j) = (b-a)/(N(j)-1);
            
        endfor
        fitting = polyfit (log(h), log(norma_err), 1);
        p = fitting(1);
        log_C = fitting(2);
        x_fitted = log(h);
        y_fitted = polyval (fitting, x_fitted);
        
        plot(log(h), log(norma_err), ";norma;", "marker","o");
        hold on
        plot(x_fitted, y_fitted, ";polyfit;");
        
        annotation ("textbox", [.6 .2 .3 .09], "string", ...
                 {strjoin({"p =", num2str(p)}), ...
                  strjoin({"log(C) =", num2str(log_C)})}, ...
                 "horizontalalignment", "center", "fitboxtotext", "off");
        xlabel("h");
        ylabel("norma do erro");
        hold off
        title("Analise assintotica do erro");
        print(strjoin( { "Figuras", strrep(diretorios(i, :), " ", ""), problema, "_tipo", int2str(t), ".png"}, "" ) );
        clf
    endfor
endfor
close