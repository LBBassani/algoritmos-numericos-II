 function grafico_solucao(u,x,y,n,m)
	z = zeros(m,n);
  I = 1;
	for j = 1 : m
		for i = 1 : n
			z(j,i) = u(I);
      I = I+1;
		end;
	end;
	%mesh(x,y,z);
  [xx,yy]= meshgrid (x,y);
  surf(xx,yy,z);
  
  xlabel("x")
  ylabel("y")
  zlabel("u(x,y)")
endfunction


