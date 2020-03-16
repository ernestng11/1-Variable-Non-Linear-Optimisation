%%***********************************************************
%% steepest: steepest descent method with Armijo line search.
%%
%%***********************************************************

  function  [x,d,iter] = steepest_Armijo(fun,x0,tol,printyes)

  if size(x0,1) < size(x0,2); x0 = x0'; end; 
  if (nargin <=4); printyes = 0; end 
  
  x = x0; 
  printyes = 1; 
  maxit = 100; 
%%
  fprintf('\n    iter   x(1)        x(2)       normgrad    step-len');
  fprintf('\n-----------------------------------------------------'); 
  t = 1;
  for iter = 0:maxit 
      [fx,grad]= feval(fun,x); 
      d = -grad;
      normd = norm(d); 
      if (printyes)
         fprintf('\n     %2.0f   %- 7.6f  %- 7.6f    %3.2e',iter,x(1),x(2),normd); 
      end
      if (normd < tol); break; end;       
      delta = grad'*d; 
      [t,iterstep] = Armijo(fun,x,fx,d,delta,1.5*t); 
      if (printyes)
         fprintf('     %4.3f %2.0f',t,iterstep); 
      end
      x = x + t*d;    
  end; 
  fprintf('\n');
%%***********************************************************

  function [t,k] = Armijo(fun,x0,fx0,d,delta,t0)
  
  beta = 0.7;
  sigma = 0.4;   
  maxit = 100; 
  for k=0:maxit
      t = beta^k*t0; 
      x = x0+t*d; 
      [fx]  = feval(fun,x);
	  if (fx < fx0 + (t*sigma)*delta)
	     break; 
	  end
  end 
%%***********************************************************