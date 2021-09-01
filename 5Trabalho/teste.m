%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Testes do exercicio computacional 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function teste()
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

  n = [10, 100, 1000];

  for i = 1:length(n)
    [x,u] = pvc(a,b,n(i),tipo_a,ua,sigma_a,alfa_a,beta_a,gamma_a,tipo_b,ub,sigma_b,alfa_b,beta_b,gamma_b, p, q, r);
    plot (x, u, strjoin({";pvc n = ", int2str(n(i)), ";"}), "linestyle" , "--");
    hold on
  endfor

  resp = f(x);
  plot(x, resp, ";original;", "linestyle", "-");
  hold off
  title("Primeiro teste de validação");
  print("validacao/validacao_1.png");
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

  n = [10, 100, 1000];

  for i = 1:length(n)
    [x,u] = pvc(a,b,n(i),tipo_a,ua,sigma_a,alfa_a,beta_a,gamma_a,tipo_b,ub,sigma_b,alfa_b,beta_b,gamma_b, p, q, r);
    plot (x, u, strjoin({";pvc n = ", int2str(n(i)), ";"}), "linestyle" , "--");
    hold on
  endfor

  resp = f(x);
  plot(x, resp, ";original;", "linestyle", "-");
  hold off
  title("Segundo teste de validação");
  print("validacao/validacao_2.png");
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

  n = [10, 100, 1000];

  for i = 1:length(n)
    [x,u] = pvc(a,b,n(i),tipo_a,ua,sigma_a,alfa_a,beta_a,gamma_a,tipo_b,ub,sigma_b,alfa_b,beta_b,gamma_b, p, q, r);
    plot (x, u, strjoin({";pvc n = ", int2str(n(i)), ";"}), "linestyle" , "--");
    hold on
  endfor

  resp = f(x);
  plot(x, resp, ";original;", "linestyle", "-");
  hold off
  title("Terceiro teste de validação");
  print("validacao/validacao_3.png");
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

  n = [10, 100, 1000];

  for i = 1:length(n)
    [x,u] = pvc(a,b,n(i),tipo_a,ua,sigma_a,alfa_a,beta_a,gamma_a,tipo_b,ub,sigma_b,alfa_b,beta_b,gamma_b, p, q, r);
    plot (x, u, strjoin({";pvc n = ", int2str(n(i)), ";"}), "linestyle" , "--");
    hold on
  endfor

  resp = f(x);
  plot(x, resp, ";original;", "linestyle", "-");
  hold off
  title("Quarto teste de validação");
  print("validacao/validacao_4.png");
  close

  %%%%% Assinatura: function [x,u] = pvc(a,b,n,tipo_a,ua,sigma_a,alfa_a,beta_a,gamma_a,tipo_b,ub,sigma_b,alfa_b,beta_b,gamma_b);
endfunction