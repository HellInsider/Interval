## -*- texinfo -*-
## @documentencoding UTF-8
## @defmethod {@@kinterval} columns (@var{A})
##
## Return the number of columns of @var{A}.
## @seealso{@@kinterval/numel, @@kinterval/size, @@kinterval/length, @@kinterval/rows, @@kinterval/end}
## @end defmethod

function result = columns (a)

  if (nargin ~= 1)
    print_usage ();
    return
  endif

  result = columns (a.inf);

endfunction

%!assert (columns (kinterval (zeros (3, 4))), 4);
%!assert (columns (kinterval (zeros (0, 4))), 4);
%!assert (columns (kinterval (zeros (3, 0))), 0);
%!assert (columns (kinterval (zeros (3, 1))), 1);
