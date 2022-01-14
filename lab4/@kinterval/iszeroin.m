function t = iszeroin(x)
  if (not (isa (x, "kinterval")))
    x = kinterval (x);
  endif
  t = in(0,x);
endfunction

