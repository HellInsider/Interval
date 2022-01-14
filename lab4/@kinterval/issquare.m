## -*- texinfo -*-
## @documentencoding UTF-8
## @defmethod {@@kinterval} issquare (@var{A})
##
## Return true if @var{A} is a square interval matrix.
##
## @seealso{@@kinterval/isvector, @@kinterval/ismatrix, @@kinterval/isscalar}
## @end defmethod

## FIXME This function is only required, because of a regression in octave 4.0
function result = issquare (A)

  if (nargin ~= 1)
    print_usage ();
    return
  endif

  result = issquare (A.inf);

endfunction

%!assert (issquare (kinterval ([])));
%!assert (issquare (kinterval (0)));
%!assert (not (issquare (kinterval (zeros (1, 2)))));
%!assert (not (issquare (kinterval (zeros (2, 1)))));
%!assert (issquare (kinterval (zeros (5))));
