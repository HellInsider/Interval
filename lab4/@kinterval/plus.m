function z = plus (x, y)
  if (not (isa (x, "kinterval")))
    x = kinterval (x);
  endif
  if (not (isa (y, "kinterval")))
    y = kinterval (y);
  endif
  if ~all(size(x.inf) == size(y.inf))
    error("wrong dimensions");
  endif

  zinf = x.inf + y.inf;
  zsup = x.sup + y.sup;

  z = kinterval (zinf, zsup);
endfunction
