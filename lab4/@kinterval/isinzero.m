function t = isinzero(x)
  if (not (isa (x, "kinterval")))
    x = kinterval (x);
  endif
  t = iszeroin(dual(x));
endfunction