## -*- texinfo -*-
## @documentencoding UTF-8
## @defmethod {@@kinterval} isrow (@var{A})
##
## Return true if @var{A} is an interval row vector.
##
## @seealso{@@kinterval/isvector, @@kinterval/iscolumn}
## @end defmethod

## FIXME This function is only required, because of a regression in octave 4.0
function result = isrow (A)

  if (nargin ~= 1)
    print_usage ();
    return
  endif

  result = isrow (A.inf);

endfunction

%!assert (not (isrow (kinterval ([]))));
%!assert (isrow (kinterval (0)));
%!assert (isrow (kinterval (zeros (1, 2))));
%!assert (not (isrow (kinterval (zeros (2, 1)))));
%!assert (not (isrow (kinterval (zeros (5)))));
