function result = mid(x)
  if (not (isa (x, "kinterval")))
    x = kinterval (x);
  endif

  z = x ./ 2;
  result = z.inf + z.sup;
endfunction