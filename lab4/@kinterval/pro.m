function z = pro(x)
  if (not (isa (x, "kinterval")))
    x = kinterval (x);
  endif
  zinf = min(x.inf, x.sup);
  zsup = max(x.inf, x.sup);

  z = kinterval(zinf, zsup);
endfunction