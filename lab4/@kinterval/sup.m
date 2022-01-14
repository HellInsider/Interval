function result = sup (x)

  if (nargin ~= 1)
    print_usage ();
    return
  endif

  result = x.sup;

endfunction

