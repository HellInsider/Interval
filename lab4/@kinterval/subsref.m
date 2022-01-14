function A = subsref (A, S)
  #printf("subsref\n");

  if (nargin ~= 2)
    print_usage ();
    return
  endif

  switch S(1).type
    case "()"
      A.inf = subsref (A.inf, S(1));
      A.sup = subsref (A.sup, S(1));
    case "{}"
      error ("kinterval cannot be indexed with {}")
    case "."
      if (not (any (strcmp (S(1).subs, methods (A)))))
        error (["kinterval property ‘", S(1).subs, "’ is unknown"])
      endif
      functionname = ["@kinterval", filesep(), S(1).subs];
      if (nargin (functionname) ~= 1)
        error (["‘", S(1).subs, "’ is not a valid kinterval property"])
      endif
      A = feval (S(1).subs, A);
    otherwise
      error ("invalid subscript type")
  endswitch

  if (numel (S) > 1)
    A = subsref (A, S(2 : end));
  endif

endfunction
