function z = and(x, y)
  if (not (isa (x, "kinterval")))
    x = kinterval (x);
  endif
  if (not (isa (y, "kinterval")))
    y = kinterval (y);
  endif
  if ~all(size(x.inf) == size(y.inf))
    error("wrong dimensions");
  endif
  z = kinterval (max(x.inf, y.inf), min(x.sup, y.sup));
endfunction