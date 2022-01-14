## -*- texinfo -*-
## @documentencoding UTF-8
## @defop Method {@@kinterval} ctranspose (@var{X})
## @defopx Operator {@@kinterval} {@var{X}'}
##
## Return the complex conjugate transpose of @var{X}.
##
## Since intervals are only defined for real numbers, this function and
## @code{transpose} are equivalent.
##
## @seealso{@@kinterval/transpose}
## @end defop

function result = ctranspose (x)

  if (nargin ~= 1)
    print_usage ();
    return
  endif

  result = transpose (x);

endfunction

