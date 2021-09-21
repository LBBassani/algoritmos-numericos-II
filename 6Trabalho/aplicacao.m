function [bx,by,gamma,fun,kappa,left,gleft,right,gright,bottom,gbottom,top,gtop] = aplicacao (x,y,n,m, app)
  switch app
    case 1
      [bx,by,gamma,fun,kappa,left,gleft,right,gright,bottom,gbottom,top,gtop] = aplicacao1 (x,y,n,m);
    case 2
      [bx,by,gamma,fun,kappa,left,gleft,right,gright,bottom,gbottom,top,gtop] = aplicacao2 (x,y,n,m);
    case 3
      [bx,by,gamma,fun,kappa,left,gleft,right,gright,bottom,gbottom,top,gtop] = aplicacao3 (x,y,n,m,1);
    case 4
      [bx,by,gamma,fun,kappa,left,gleft,right,gright,bottom,gbottom,top,gtop] = aplicacao3 (x,y,n,m,2);
  endswitch
endfunction


function [bx,by,gamma,fun,kappa,left,gleft,right,gright,bottom,gbottom,top,gtop] = aplicacao1 (x,y,n,m)
 I = 1;
 for i = 1:n
     for j = 1:m
         bx(I) = 0.0;
         by(I) = 0.0;
         gamma(I) = 0.0;
         fun(I) = 0.0;
         I++;
     endfor
 endfor
 kappa = 1.0;

 left = 1;
 gleft = 10.0;

 right = 1;
 gright = 10.0;

 bottom = 1; 
 gbottom = 10.0;

 top = 1;
 gtop = 10.0;

endfunction 

function [bx,by,gamma,fun,kappa,left,gleft,right,gright,bottom,gbottom,top,gtop] = aplicacao2 (x,y,n,m)
 I = 1;
 for i = 1:n
     for j = 1:m
         bx(I) = 0.0;
         by(I) = 0.0;
         gamma(I) = 0.0;
         fun(I) = 0.0;
         I++;
     endfor
 endfor
 kappa = 1.0;

 left = 1;
 gleft = 10.0;

 right = 1;
 gright = 10.0;

 bottom = 1; 
 gbottom = 100.0;

 top = 1;
 gtop = 100.0;

endfunction 

function [bx,by,gamma,fun,kappa,left,gleft,right,gright,bottom,gbottom,top,gtop] = aplicacao3 (x,y,n,m,ver)
  kappa = 1.0;
  
  u = @(x, y) 10 .* (x .* y) .* (1 - x) .* (1 - y) .* e.^(x.^4.5);
  u_x = @(x, y) -45 .* e.^(x.^4.5) .* (0.222222 - 0.444444 .* x + x.^4.5 - x.^5.5) .* (-1 + y) .* y;
  u_x_x = @(x, y) 337.5 .* e.^(x.^4.5) .* (0.0592592 - 0.733333 .* x.^3.5 + x.^4.5 - 0.6 .* x.^8 + 0.6 .* x.^9) .* (-1 + y) .* y;
  u_y = @(x, y) 10 .* e.^(x.^4.5) .* (-1 + x) .* x .* (-1 + 2 .* y);
  u_y_y = @(x, y) 20 .* e.^(x.^4.5) .* (-1 + x) .* x;
  
  [X, Y] = meshgrid(x, y);
  ureal = u(X, Y);
  
 I = 1;
 for j = 1:m
     for i = 1:n
         bx(I) = 1.0;
         by(I) = 20.0*y(j);
         gamma(I) = 1.0;
         I++;
     endfor
 endfor
 
 fun = -kappa * (u_x_x(X, Y)'(:) + u_y_y(X, Y)'(:)) + bx .* u_x(X, Y)'(:) + by .* u_y(X, Y)'(:) + ureal'(:);
 
 switch ver
   case 1
       
     left = 1;
     gleft = 0.0;

     right = 1;
     gright = 0.0;

     bottom = 1; 
     gbottom = 0.0;

     top = 1;
     gtop = 0.0;
   
    case 2
   
     x_index = @(I) rem(I, n);
     y_index = @(I) floor(I/n) + 1;
          

     I = [1:n:(m-1)*n+1];
     gleft = sparse(n*m,1);
          
     x_i = x(x_index(I));
     y_i = y(y_index(I));
     du = u_x(x_i, y_i);
     
     left = 2;
     gleft(I) = du*kappa;
   
     right = 1;
     gright = 0.0;

     bottom = 1; 
     gbottom = 0.0;

     top = 1;
     gtop = 0.0;

 endswitch 
endfunction


##
##function [bx,by,gamma,fun,kappa,left,gleft,right,gright,bottom,gbottom,top,gtop] = aplicacao (x,y,n,m)
## I = 1;
## for i = 1:n
##     for j = 1:m
##         bx(I) = ...;
##         by(I) = ...;
##         gamma(I) = ...;
##         fun(I) = ...;
##         I++;
##     endfor
## endfor
## kappa = ...;
##
## left = ... ;
## gleft = ...;
##
## right = ...;
## gright = ...;
##
## bottom = ...; 
## gbottom = ...;
##
## top = ...;
## gtop = ...;
##
##endfunction 
##
