function result = isscalar (A)

  if (nargin ~= 1)
    print_usage ();
    return
  endif

  result = isscalar (A.inf);

endfunction
