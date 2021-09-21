 function [a,b,c,d,e] = coeficientes (hx,hy,kappa,bx,by,gamma,n,m);

   N = n*m;
   
   for I=1:N
      a(I) = gamma(I) + 2.0*kappa*(1.0/(hx^2) + 1.0/(hy^2));
      b(I) = - kappa/(hx^2) - (bx(I)/(2.0*hx));
      c(I) = - kappa/(hx^2) + (bx(I)/(2.0*hx));
      d(I) = - kappa/(hy^2) - (by(I)/(2.0*hy));
      e(I) = - kappa/(hy^2) + (by(I)/(2.0*hy));
   endfor
 endfunction
