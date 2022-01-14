function z = minus (x, y)
  if (not (isa (x, "kinterval")))
    x = kinterval (x);
  endif

  if (not (isa (y, "kinterval")))
    y = kinterval (y);
  endif

  if ~all(size(x.inf) == size(y.inf))
    error("wrong dimensions");
  endif

  zinf = x.inf - y.sup;
  zsup = x.sup - y.inf;

  z = kinterval (zinf, zsup);
endfunction
