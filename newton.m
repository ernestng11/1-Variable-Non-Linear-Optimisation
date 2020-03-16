%%************************************************************
%% newton: 
%%
%% [x,iter] = newton(fun,x0,tol,funparms,printyes); 
%% fun = function file where df(x) and d2f(x) are evaluated;
%% x0 = initial iterate;
%% funparms = parameters used in fun
%% printyes = 1 if want to print information at each iterations.
%%
%% Example: [x,iter] = newton('newtonfun',0.4,1e-6,[],1); 
%%************************************************************

  function  [x,iter] = newton(fun,x0,tol,funparms,printyes); 

   if (nargin < 5); printyes = 0; end
   if (nargin < 4); funparms = []; end
   
   if (printyes)
      fprintf('\n iter     xk        f(xk)       df(xk)\n')
      fprintf('----------------------------------------\n')
   end
   x = x0; 
   maxit = 50; 
   for iter = 0:maxit
       [f,df,d2f] = feval(fun,x,funparms);    
       if (printyes) 
          fprintf('  %2.0f    %- 5.4f   %- 5.4f    %- 3.2e\n',iter,x,f,df);
       end
       if (abs(df) < tol); break; end; 
       x = x - df/d2f;        
   end
%%************************************************************



