## -*- texinfo -*-
## @documentencoding UTF-8
## @defmethod {@@kinterval} ismatrix (@var{A})
##
## Return true if @var{A} is a 2-D interval array.
##
## @seealso{@@kinterval/isscalar, @@kinterval/isvector}
## @end defmethod

## Author: Oliver Heimlich
## Keywords: interval
## Created: 2014-10-29

## FIXME: One happy day this function can be removed, because bug #42422 has
## been solved with GNU Octave 4.0.
function result = ismatrix (A)

  if (nargin ~= 1)
    print_usage ();
    return
  endif

  result = ismatrix (A.inf);

endfunction

%!assert (ismatrix (kinterval ([])));
%!assert (ismatrix (kinterval (0)));
%!assert (ismatrix (kinterval (zeros (3, 1))));
%!assert (ismatrix (kinterval (zeros (1, 4))));
%!assert (ismatrix (kinterval (zeros (3, 4))));
