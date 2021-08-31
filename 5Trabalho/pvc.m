function [x,u]=pvc(a,b,n,tipo_a,ua,sigma_a,alfa_a,beta_a,gamma_a,tipo_b,ub,sigma_b,alfa_b,beta_b,gamma_b, p, q, r);

    % define variaveis principais
    h = (b-a)/(n-1);
    x = linspace(a,b,n);
    x = x';
    f = zeros(n,1);
    u = zeros(n,1);
    A = zeros(n,n);

    % monta sistema Au=f
    for i = 1:n
        y = x(i); % Valor para aplicar as funções
        
        if i == 1
          b1 = (1/realpow(h, 2)) - p(y)/(2*h); % b_1 : primeiro b não entra na matriz
        else
          A(i,i-1) = (1/realpow(h, 2)) - p(y)/(2*h); % b_i
        endif
        
        A(i,i) = q(y) - (2/realpow(h, 2));  % a_i
        
        if i == n
          cn       = (1/realpow(h, 2)) + p(y)/(2*h); % c_n : ultimo c não entra na matriz
        else
          A(i,i+1) = (1/realpow(h, 2)) + p(y)/(2*h); % c_i
        endif
        
        f(i) = r(y); % r_i
    endfor

    % aplica condicoes de contorno em x=a;
    switch tipo_a
        case 1  % valor prescrito
            A(1, 1) = 1;
            A(1, 2) = 0;
            f(1) = ua;
        case 2 % derivada prescrita
            A(1,1) = A(1,1) + b1;
            f(1) = f(1) + b1*h*sigma_a;
        case 3 % misto
            A(1,1) = A(1,1) + b1* (1 + (h*beta_a)/alfa_a);
            f(1) = f(1) + (b1*h*gamma_a)/alfa_a;
        otherwise
               printf("Erro na Condicao de contorno"); 
    end

    % aplica condicoes de contorno em x=b;
    switch tipo_b
        case 1 % valor prescrito
            A(n, n - 1) = 0;
            A(n, n) = 1;
            f(n) = ub;
        case 2 % derivada prescrita
            A(n, n - 1) = A(n,n) + cn;
            f(n) = f(n) - cn*h*sigma_b;
        case 3 % misto
            A(n,n) = A(n,n)  + cn* (1 - (h*beta_b)/alfa_b);
            f(n) = f(n) - (cn*h*gamma_b)/alfa_b;
        otherwise
            printf("Erro na Condicao de contorno"); 
    end 

    % Resolve sistema linear
    u = A\f;
      
endfunction       
       
