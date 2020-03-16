%%************************************************************
%% bisection: 
%%
%% [x,iter] = bisection(fun,a,b,tol,funparms,printyes); 
%% fun = function file where df(x) is evaluated;
%% a,b = initial interval;
%% funparms = parameters used in fun
%% printyes = 1 if want to print information at each iterations.
%%
%%  [x,iter,flag] = bisection('bisectionfun',-1,1,0.01,[],1);
%%************************************************************

  function  [x,iter,flag] = bisection(fun,a,b,tol,funparms,printyes); 

   if (nargin <= 5); printyes = 0; end
   if (nargin <= 4); funparms = []; end
   
   flag = 0; 
   maxit = ceil(log((b-a)/(tol+eps))/log(2));
   dfa = feval(fun,a,funparms);
   dfb = feval(fun,b,funparms);
   if (sign(dfa)*sign(dfb) > 0)
      iter = 0; flag = 1;
      x = (a+b)/2;
      fprintf(' dfa, dfb have same sign');  return;
   end
   if (printyes)
      fprintf('   k      ak       bk        xk       df(xk)\n') 
      fprintf('---------------------------------------------\n')
   end
%%
   for iter = 1:maxit
       x = (a+b)/2; 
       dfx = feval(fun,x,funparms);
       if (printyes)
           fprintf('  %2.0f   %- 5.4f   %- 5.4f   %- 5.4f   %- 5.4f\n',...
        	     iter,a,b,x,dfx);
       end
       if (sign(dfx)*sign(dfb) <= 0)
          a = x; dfa = dfx; 
       elseif (sign(dfx)*sign(dfa) <= 0)
          b = x; dfb = dfx; 
       end
       if (b-a) < tol;  break; end;    
   end 
   if (printyes)
      fprintf('---------------------------------------------\n')
   end
%%************************************************************
