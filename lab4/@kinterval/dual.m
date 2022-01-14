function z = dual(x)
  if (not (isa (x, "kinterval")))
    x = kinterval (x);
  endif
  z = kinterval(x.sup, x.inf);
endfunction