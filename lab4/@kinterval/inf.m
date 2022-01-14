function result = inf (x)

  if (nargin ~= 1)
    print_usage ();
    return
  endif

  result = x.inf;

endfunction
