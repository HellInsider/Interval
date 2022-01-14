## -*- texinfo -*-
## @documentencoding UTF-8
## @defmethod {@@kinterval} reshape (@var{A}, @var{M}, @var{N}, ...)
## @defmethodx {@@kinterval} reshape (@var{X}, [@var{M} @var{N}, ...])
## @defmethodx {@@kinterval} reshape (@var{X}, ..., @var{[]}, ...)
##
## Return an interval matrix with the specified dimensions (M, N, ...) whose
## elements are taken from the interval matrix @var{A}.  The elements of the
## matrix are accessed in column-major order (like Fortran arrays are stored).
##
## Note that the total number of elements in the original matrix
## (@code{prod (size (@var{A}))}) must match the total number of elements in
## the new matrix (@code{prod ([@var{M} @var{N}])}).
##
## A single dimension of the return matrix may be left unspecified and
## Octave will determine its size automatically.  An empty matrix ([])
## is used to flag the unspecified dimension.
##
## @example
## @group
## reshape (kinterval (1 : 6), 2, 3)
##   @result{} ans = 2×3 interval matrix
##          [1]   [3]   [5]
##          [2]   [4]   [6]
## @end group
## @end example
## @seealso{@@kinterval/resize, @@kinterval/cat, @@kinterval/postpad, @@kinterval/prepad}
## @end defmethod

function a = reshape (a, varargin)

  if (not (isa (a, "kinterval")))
    print_usage ();
    return
  endif

  a.inf = reshape (a.inf, varargin{:});
  a.sup = reshape (a.sup, varargin{:});

endfunction

%!assert (reshape (kinterval (1 : 6), 2, 3) == kinterval (reshape (1 : 6, 2, 3)));
%!assert (reshape (kinterval (1 : 24), 2, [], 4) == kinterval (reshape (1 : 24, 2, 3, 4)));
