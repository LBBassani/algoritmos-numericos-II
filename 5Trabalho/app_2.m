%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Resfriador unidimensional %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function app_2()
    % Valores do problema
    u_0 = 160;
    u_ref = 70;
    K = 0.001;
    T = 0.1;
    W = 10;
    L = 1;
    c_ref = [0.1; 0.01; 0.001; 0.0001];
    n = [10; 50; 100];
    
    for j = 1:length(c_ref)    
        
        C = ((2*W + 2*T)/T*W) * c_ref(j);
        r = @(x) -(C*u_ref)/K;
        
        % y'' + Cy = (-1/K)f(x)
        p = @(x) 0;
        q = @(x) -C/K;
        
        tipo_a = 1;
        a = 0;
        ua = u_0;
        sigma_a = alfa_a = beta_a = gamma_a = 0;
        
        tipo_b = 3;
        b = L;
        alfa_b = K;
        beta_b = c_ref(j);
        gamma_b = c_ref(j)*u_ref;
        ub = sigma_b = 0;
        
        for i = 1:length(n)
            [x,u] = pvc(a,b,n(i),tipo_a,ua,sigma_a,alfa_a,beta_a,gamma_a,tipo_b,ub,sigma_b,alfa_b,beta_b,gamma_b, p, q, r);
            figure 1
            title(strjoin({"Gráfico da solução para c_{ref} = ", num2str(c_ref(j), "%.1e")}))
            plot (x, u, strjoin({";pvc n = ", int2str(n(i)), ";"}), "linestyle" , "--");
            hold on
        endfor
        print(strjoin({"app_2/a_", num2str(log10(c_ref(j))), ".png"}, ""))
        hold off
        figure 2
        title("Gráfico da solução para n = 1000")
        [x,u] = [x,u] = pvc(a,b,1000,tipo_a,ua,sigma_a,alfa_a,beta_a,gamma_a,tipo_b,ub,sigma_b,alfa_b,beta_b,gamma_b, p, q, r);
        plot (x, u, strjoin({";pvc c_{ref} = ", num2str(c_ref(j), "%.1e") , ";"}), "linestyle" , "--");
        hold on
     endfor
     
     print("app_2/a.png")
endfunction