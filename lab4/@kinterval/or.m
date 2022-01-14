    function z = or(x, y)
      if (not (isa (x, "kinterval")))
        x = kinterval (x);
      endif
      if (not (isa (y, "kinterval")))
        y = kinterval (y);
      endif
      if ~all(size(x.inf) == size(y.inf))
        error("wrong dimensions");
      endif
      z = kinterval (min(x.inf, y.inf), max(x.sup, y.sup));
    endfunction