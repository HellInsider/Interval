## -*- texinfo -*-
## @documentencoding UTF-8
## @defop Method {@@infsup} transpose (@var{X})
## @defopx Operator {@@infsup} {@var{X}.'}
##
## Return the transpose of interval matrix or vector @var{X}.
##
## @example
## @group
## infsup (zeros (1, 3), ones (1, 3)) .'
##   @result{} ans = 3×1 interval vector
##      [0, 1]
##      [0, 1]
##      [0, 1]
## @end group
## @end example
## @seealso{@@infsup/ctranspose}
## @end defop

function x = transpose (x)

  if (nargin ~= 1)
    print_usage ();
    return
  endif

  x.inf = transpose (x.inf);
  x.sup = transpose (x.sup);

endfunction

