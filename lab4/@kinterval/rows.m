## -*- texinfo -*-
## @documentencoding UTF-8
## @defmethod {@@kinterval} rows (@var{A})
##
## Return the number of rows of @var{A}.
## @seealso{@@kinterval/numel, @@kinterval/size, @@kinterval/length, @@kinterval/columns, @@kinterval/end}
## @end defmethod

function result = rows (a)

  if (nargin ~= 1)
    print_usage ();
    return
  endif

  result = rows (a.inf);

endfunction

%!assert (rows (kinterval (zeros (3, 4))), 3);
