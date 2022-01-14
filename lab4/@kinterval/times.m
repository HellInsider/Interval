   function z = times(x, y)
     
      if (not (isa (x, "kinterval")))
        x = kinterval (x);
      endif
      if (not (isa (y, "kinterval")))
        y = kinterval (y);
      endif
      
      if isscalar (x) && ~isscalar (y) 
          x.inf = x.inf * ones(size(y.inf));
          x.sup = x.sup * ones(size(y.inf));
      end

      if ~isscalar (x) && isscalar (y) 
          y.inf = y.inf * ones(size(x.inf));
          y.sup = y.sup * ones(size(x.inf));
      end
     
      if ~all(size(x.inf) == size(y.inf))
        error("wrong dimensions");
      endif
    
      xinfm = max(-x.inf, 0);
      xinfp = max(x.inf, 0);
      yinfm = max(-y.inf, 0);
      yinfp = max(y.inf, 0);

      xsupm = max(-x.sup, 0);
      xsupp = max(x.sup, 0);
      ysupm = max(-y.sup, 0);
      ysupp = max(y.sup, 0);
      
      zinf = max(xinfp .* yinfp, xsupm .* ysupm) - max(xsupp .* yinfm, xinfm .* ysupp);
      zsup = max(xsupp .* ysupp, xinfm .* yinfm) - max(xinfp .* ysupm, xsupm .* yinfp);
      
      z = kinterval(zinf, zsup);
     
 endfunction 
