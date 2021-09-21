function [bx,by,gamma,fun,kappa,left,gleft,right,gright,bottom,gbottom,top,gtop] = aplicacao (x,y,n,m, app)
  switch app
    case 1
      [bx,by,gamma,fun,kappa,left,gleft,right,gright,bottom,gbottom,top,gtop] = aplicacao1 (x,y,n,m);
    case 2
      [bx,by,gamma,fun,kappa,left,gleft,right,gright,bottom,gbottom,top,gtop] = aplicacao2 (x,y,n,m);
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
