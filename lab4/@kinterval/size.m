function varargout = size (a, dim)

  if (nargin == 0 || nargin > 2)
    print_usage ();
    return
  endif

  if (nargin == 1)
     varargout = cell (1, max (1, nargout));
     [varargout{:}] = size (a.inf);
  else
     if (nargout > 1)
       print_usage ();
       return
     endif
     varargout{1} = size (a.inf, dim);
  endif

endfunction

