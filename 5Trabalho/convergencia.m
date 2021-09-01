function convergencia()
    
    % Análise de convergencia para a norma do máximo considerando os problemas com solução conhecida
    n = [4, 8, 16, 32];
    
    
    % Primeiro teste: 
    % y'' + xy' + y = 4x^3 + 4x + 1
    % p(x) = x, q(x) = 1, r(x) = 4x^3 + 4x + 1, y(x) = x^3 - x + 1
    p = @(x) x;
    q = @(x) 1;
    r = @(x) 4*realpow(x,3) + 4*x + 1;
    f = @(x) realpow(x,3) - x + 1;

    % y(0) = 1 {Condição tipo 1}; ua = 1
    tipo_a = 1;
    a = 0;
    ua = 1;
    sigma_a = alfa_a = beta_a = gamma_a = 0;

    % 2y'(1) + y(1) = 5 {Condição tipo 3}; alpha_b = 2, beta_b = 1, gama_b = 5
    tipo_b = 3;
    b = 1;
    alfa_b = 2;
    beta_b = 1;
    gamma_b = 5;
    ub = sigma_b = 0;
    
    norma_err = zeros(1, length(n));
    h = zeros(1, length(n));

    for i = 1:length(n)
        [x,u] = pvc(a,b,n(i),tipo_a,ua,sigma_a,alfa_a,beta_a,gamma_a,tipo_b,ub,sigma_b,alfa_b,beta_b,gamma_b, p, q, r);
        resp = f(x);
        err = resp - u;
        norma_err(i) = norm(err, Inf);
        h(i) = (b-a)/(n(i)-1);
    endfor
    fitting = polyfit (log(h), log(norma_err), 1);
    p = fitting(1);
    log_C = fitting(2);
    x = log(h);
    y = polyval (fitting, x);
    
    title("Primeiro teste de validação");
    plot(log(h), log(norma_err), ";norma;", "marker","o");
    hold on
    plot(x, y, ";polyfit;");
    
    annotation ("textbox", [.6 .2 .3 .09], "string", ...
             {strjoin({"p =", num2str(p)}), ...
              strjoin({"log(C) =", num2str(log_C)})}, ...
             "horizontalalignment", "center", "fitboxtotext", "off");
    xlabel("h");
    ylabel("norma do erro");
    hold off
    print("validacao/convergencia_1.png");
    clf
    
    % Segundo teste: 
    % y'' = 1
    % p(x) = 0, q(x) = 0, r(x) = 1, y(x) = x^2/2 - x + 1
    p = @(x) 0;
    q = @(x) 0;
    r = @(x) 1;
    f = @(x) realpow(x, 2)/2 - x + 1;

    % y(0) = 1 {Condição tipo 1}; ua = 1
    tipo_a = 1;
    a = 0;
    ua = 1;
    sigma_a = alfa_a = beta_a = gamma_a = 0;

    % y(1) = 1/2 {Condição tipo 1}; ub = 0
    tipo_b = 1;
    b = 1;
    ub = 1/2;
    sigma_b = alfa_b = beta_b = gamma_b = 0;
    
    norma_err = zeros(1, length(n));
    h = zeros(1, length(n));

    for i = 1:length(n)
        [x,u] = pvc(a,b,n(i),tipo_a,ua,sigma_a,alfa_a,beta_a,gamma_a,tipo_b,ub,sigma_b,alfa_b,beta_b,gamma_b, p, q, r);
        resp = f(x);
        err = resp - u;
        norma_err(i) = norm(err, Inf);
        h(i) = (b-a)/(n(i)-1);
    endfor
    fitting = polyfit (log(h), log(norma_err), 1);
    p = fitting(1);
    log_C = fitting(2);
    x = log(h);
    y = polyval (fitting, x);
    
    title("Segundo teste de validação");
    plot(log(h), log(norma_err), ";norma;", "marker","o");
    hold on
    plot(x, y, ";polyfit;");
    
    annotation ("textbox", [.6 .2 .3 .09], "string", ...
             {strjoin({"p =", num2str(p)}), ...
              strjoin({"log(C) =", num2str(log_C)})}, ...
             "horizontalalignment", "center", "fitboxtotext", "off");
    xlabel("h");
    ylabel("norma do erro");
    hold off
    print("validacao/convergencia_2.png");
    clf
    
    % Terceiro teste: 
    % y'' = x + 1
    % p(x) = 0, q(x) = 0, r(x) = x + 1, u(x) = x^3/6 + x^2/2 - x/6 - 1/2
    p = @(x) 0;
    q = @(x) 0;
    r = @(x) 1;
    f = @(x) realpow(x, 3)/6 + realpow(x, 2)/2 - x/6 - 1/2;

    % y(-1) = 0 {Condição tipo 1}; ua = 0
    tipo_a = 1;
    a = -1;
    ua = 0;
    sigma_a = alfa_a = beta_a = gamma_a = 0;

    % y(1) = 0 {Condição tipo 1}; ub = 0
    tipo_b = 1;
    b = 1;
    ub = 0;
    sigma_b = alfa_b = beta_b = gamma_b = 0;

    norma_err = zeros(1, length(n));
    h = zeros(1, length(n));

    for i = 1:length(n)
        [x,u] = pvc(a,b,n(i),tipo_a,ua,sigma_a,alfa_a,beta_a,gamma_a,tipo_b,ub,sigma_b,alfa_b,beta_b,gamma_b, p, q, r);
        resp = f(x);
        err = resp - u;
        norma_err(i) = norm(err, Inf);
        h(i) = (b-a)/(n(i)-1);
    endfor
    fitting = polyfit (log(h), log(norma_err), 1);
    p = fitting(1);
    log_C = fitting(2);
    x = log(h);
    y = polyval (fitting, x);
    
    title("Terceiro teste de validação");
    plot(log(h), log(norma_err), ";norma;", "marker","o");
    hold on
    plot(x, y, ";polyfit;");
    
    annotation ("textbox", [.6 .2 .3 .09], "string", ...
             {strjoin({"p =", num2str(p)}), ...
              strjoin({"log(C) =", num2str(log_C)})}, ...
             "horizontalalignment", "center", "fitboxtotext", "off");
    xlabel("h");
    ylabel("norma do erro");
    hold off
    print("validacao/convergencia_3.png");
    clf
   
    % Quarto teste
    % y'' - (1/2)y' + y = x^2 + 1/2
    % p(x) = - 1/2, q(x) = 1, r(x) = x^2 + 1/2
    % y(x) = x^2 + x - 1
    p = @(x) -1/2;
    q = @(x) 1;
    r = @(x) realpow(x,2) + 1/2;
    f = @(x) realpow(x,2) + x - 1;

    % Teste: u(0) = -1 e u(1) = 1
    tipo_a = 1;
    a = 0;
    ua = -1;
    sigma_a = alfa_a = beta_a = gamma_a = 0;

    tipo_b = 1;
    b = 1;
    ub = 1;
    sigma_b = alfa_b = beta_b = gamma_b = 0;

    norma_err = zeros(1, length(n));
    h = zeros(1, length(n));

    for i = 1:length(n)
        [x,u] = pvc(a,b,n(i),tipo_a,ua,sigma_a,alfa_a,beta_a,gamma_a,tipo_b,ub,sigma_b,alfa_b,beta_b,gamma_b, p, q, r);
        resp = f(x);
        err = resp - u;
        norma_err(i) = norm(err, Inf);
        h(i) = (b-a)/(n(i)-1);
    endfor
    fitting = polyfit (log(h), log(norma_err), 1);
    p = fitting(1);
    log_C = fitting(2);
    x = log(h);
    y = polyval (fitting, x);
    
    title("Quarto teste de validação");
    plot(log(h), log(norma_err), ";norma;", "marker","o");
    hold on
    plot(x, y, ";polyfit;");
    
    annotation ("textbox", [.6 .2 .3 .09], "string", ...
             {strjoin({"p =", num2str(p)}), ...
              strjoin({"log(C) =", num2str(log_C)})}, ...
             "horizontalalignment", "center", "fitboxtotext", "off");
    xlabel("h");
    ylabel("norma do erro");
    hold off
    print("validacao/convergencia_4.png");
    
    close
    
endfunction