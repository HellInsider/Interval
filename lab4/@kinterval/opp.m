    function z = opp(x)
      if (not (isa (x, "kinterval")))
        x = kinterval (x);
      endif      
	    z = kinterval(-x.inf, -x.sup);
    endfunction