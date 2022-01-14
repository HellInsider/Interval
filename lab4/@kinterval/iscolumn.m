## -*- texinfo -*-
## @documentencoding UTF-8
## @defmethod {@@kinterval} iscolumn (@var{A})
##
## Return true if @var{A} is an interval column vector.
##
## @seealso{@@kinterval/isvector, @@kinterval/isrow}
## @end defmethod

## FIXME This function is only required, because of a regression in octave 4.0
function result = iscolumn (A)

  if (nargin ~= 1)
    print_usage ();
    return
  endif

  result = iscolumn (A.inf);

endfunction

%!assert (not (iscolumn (kinterval ([]))));
%!assert (iscolumn (kinterval (0)));
%!assert (not (iscolumn (kinterval (zeros (1, 2)))));
%!assert (iscolumn (kinterval (zeros (2, 1))));
%!assert (not (iscolumn (kinterval (zeros (5)))));
