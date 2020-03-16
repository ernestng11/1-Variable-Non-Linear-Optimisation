%%***********************************************************
%% steepest: steepest descent method with exact line search.
%%
%%***********************************************************

  function  [x,d,iter] = steepest(fun,x0,tol,printyes); 

  if size(x0,1) < size(x0,2); x0 = x0'; end; 
  if (nargin <=4); printyes = 0; end 
  
  x = x0;
  tol2 = min(1e-2*tol,1e-4); 
  printyes = 1; 
  maxit = 100; 
%%
  fprintf('\n    iter   x(1)        x(2)       normgrad    step-len');
  fprintf('\n-----------------------------------------------------'); 
  for iter = 0:maxit 
      [fx,grad]= feval(fun,x); 
      d = -grad;
      normd = norm(d); 
      if (printyes)
         fprintf('\n     %2.0f   %- 6.6f  %- 6.6f    %3.2e',iter,x(1),x(2),normd); 
      end
      if (normd < tol); break; end;       
      funparms{1} = fun; funparms{2} = x; funparms{3} = d;       
      [t,dummy] = bisection('steplength',0,10,tol2,funparms,0); 
      if (printyes)
         fprintf('     %4.3f',t); 
      end
      x = x + t*d;    
  end; 
  fprintf('\n');
%%***********************************************************
