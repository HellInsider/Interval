## -*- texinfo -*-
## @documentencoding UTF-8
## @defmethod {@@kinterval} ndims (@var{A})
##
## Return the number of dimensions of @var{A}.
##
## For any array, the result will always be larger than or equal to 2.
## Trailing singleton dimensions are not counted.
## @seealso{@@kinterval/length, @@kinterval/numel, @@kinterval/size, @@kinterval/rows, @@kinterval/columns}
## @end defmethod

function result = ndims (a)

  if (nargin ~= 1)
    print_usage ();
    return
  endif

  result = ndims (a.inf);

endfunction

%!assert (ndims (kinterval ([])), 2);
%!assert (ndims (kinterval (0)), 2);
%!assert (ndims (kinterval (magic (3))), 2);
