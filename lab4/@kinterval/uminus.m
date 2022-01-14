function x = uminus (x)

  if (nargin ~= 1)
    print_usage ();
    return
  endif

  ## This also works for empty / entire intervals
  [x.inf, x.sup] = deal (-x.sup, -x.inf);

endfunction
