

  function h = steepestfun2(t,parms); 

   gradfun = parms{1}; 
   x = parms{2}; 
   d = parms{3}; 

   y = x+t*d; 
   [f,g] = feval(gradfun,y); 
   h = g'*d; 
