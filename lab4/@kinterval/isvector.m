## -*- texinfo -*-
## @documentencoding UTF-8
## @defmethod {@@kinterval} isvector (@var{A})
##
## Return true if @var{A} is an interval vector.
##
## Scalars (1x1 vectors) are subsets of the more general vector and
## @code{isvector} will return true for these objects as well.
## @seealso{@@kinterval/isscalar, @@kinterval/ismatrix}
## @end defmethod

## Author: Oliver Heimlich
## Keywords: interval
## Created: 2015-02-02

## FIXME This function is only required, because of regression #44498
function result = isvector (A)

  if (nargin ~= 1)
    print_usage ();
    return
  endif

  result = isvector (A.inf);

endfunction

%!assert (not (isvector (kinterval ([]))));
%!assert (isvector (kinterval (0)));
%!assert (isvector (kinterval (zeros (1, 2))));
%!assert (isvector (kinterval (zeros (2, 1))));
%!assert (not (isvector (kinterval (zeros (5)))));
