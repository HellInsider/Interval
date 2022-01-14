function t = in(x, y)
      if (not (isa (x, "kinterval")))
        x = kinterval (x);
      endif
      if (not (isa (y, "kinterval")))
        y = kinterval (y);
      endif
      if isscalar(x) && ~isscalar(y)
         x = x .* ones(size(y));
      endif
      if ~isscalar(x) && isscalar(y)
         y = y .* ones(size(x));
      endif
      if ~all(size(x.inf) == size(y.inf))
        error("wrong dimensions");
      endif
      t = x.inf >= y.inf & x.sup <= y.sup;
endfunction
