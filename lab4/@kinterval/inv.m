function z = inv(x)
  if (not (isa (x, "kinterval")))
    x = kinterval (x);
  endif

  if iszeroin(pro(x))
    error ('wrong division, pro(x) contains zero');
  else
    z = kinterval(1/x.sup, 1/x.inf);
  endif
endfunction
