%%%%%%%%%%%%%%%%%%%%%%%%% Conservação de Calor em uma haste longa e fina %%%%%%%%%%%%%%%%%%%%%%%%%%
function app_1()
  % (A)
  % Valores do problema
  K = 0.01; % m²
  T_ar = 20; % ºC
  a = 0; % m
  b = 10; % m
  T_a = 40; % ºC
  T_b = 200; % ºC

  % y'' - y' = K*T_ar
  % p(x) = -1, q(x) = 0, r(x) = K*T_ar
  p = @(x) -1;
  q = @(x) 0;
  r = @(x) K*T_ar;

  tipo_a = 1;
  ua = T_a;
  sigma_a = alfa_a = beta_a = gamma_a = 0;

  tipo_b = 1;
  ub = T_b;
  sigma_b = alfa_b = beta_b = gamma_b = 0;

  n = [10, 50, 100, 1000];

  for i = 1:length(n)
    [x,u] = pvc(a,b,n(i),tipo_a,ua,sigma_a,alfa_a,beta_a,gamma_a,tipo_b,ub,sigma_b,alfa_b,beta_b,gamma_b, p, q, r);
    figure 1
    title(strjoin({"Gráfico da solução para n = ", int2str(n(i))}))
    plot (x, u, strjoin({";pvc n = ", int2str(n(i)), ";"}), "linestyle" , "--");
    print(strjoin({"app_1/a_", int2str(n(i)), ".png"}, ""))
    figure 2
    plot (x, u, strjoin({";pvc n = ", int2str(n(i)), ";"}), "linestyle" , "--");
    hold on
  endfor

  print(strjoin({"app_1/a", ".png"}, ""))
  hold off

  % (B)
  tipo_b = 2;
  sigma_b = 0;

  n = [10, 50, 100, 1000];

  for i = 1:length(n)
    [x,u] = pvc(a,b,n(i),tipo_a,ua,sigma_a,alfa_a,beta_a,gamma_a,tipo_b,ub,sigma_b,alfa_b,beta_b,gamma_b, p, q, r);
    figure 1
    title(strjoin({"Gráfico da solução para n = ", int2str(n(i))}))
    plot (x, u, strjoin({";pvc n = ", int2str(n(i)), ";"}), "linestyle" , "--");
    print(strjoin({"app_1/b_", int2str(n(i)), ".png"}, ""))
    figure 2
    plot (x, u, strjoin({";pvc n = ", int2str(n(i)), ";"}), "linestyle" , "--");
    hold on
  endfor

  print(strjoin({"app_1/b", ".png"}, ""))
  hold off
endfunction