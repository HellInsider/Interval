function result = sum (x, dim)
    
  if (nargin > 2)
    print_usage ();
    return
  endif
  
  if (nargin < 2)
    ## Try to find non-singleton dimension
    dim = find (size (x.inf) ~= 1, 1);
    if (isempty (dim))
      dim = 1;
    endif
  endif

  ## Short-circuit
  if (size (x.inf, dim) == 1)
    result = x;
    return
  endif
  
  zinf = sum (x.inf, dim);
  zsup = sum(x.sup, dim);

  result = kinterval(zinf, zsup);

endfunction
