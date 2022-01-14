function z = abs(x)
  if (not (isa (x, "kinterval")))
    x = kinterval (x);
  endif
  z = max(abs(x.inf), abs(x.sup));
endfunction