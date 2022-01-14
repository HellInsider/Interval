    function t = ispositive(x)
      if (not (isa (x, "kinterval")))
        x = kinterval (x);
      endif
    	t = x.inf >= 0 & x.sup >= 0;
    endfunction