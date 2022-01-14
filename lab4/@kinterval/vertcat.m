## -*- texinfo -*-
## @documentencoding UTF-8
## @defop Method {@@kinterval} vertcat (@var{ARRAY1}, @var{ARRAY2}, @dots{})
## @defopx Operator {@@kinterval} {[@var{ARRAY1}; @var{ARRAY2}; @dots{}]}
##
## Return the vertical concatenation of interval array objects along
## dimension 1.
##
## @example
## @group
## a = kinterval (2, 5);
## [a; a; a]
##   @result{} ans = 3×1 interval vector
##      [2, 5]
##      [2, 5]
##      [2, 5]
## @end group
## @end example
## @seealso{@@kinterval/horzcat}
## @end defop

function result = vertcat (varargin)

  result = cat (1, varargin{:});

endfunction

%!assert (vertcat (kinterval (1), kinterval (2)) == kinterval (vertcat (1, 2)));
%!test
%! # from the documentation string
%! a = kinterval (2, 5);
%! assert (vertcat (a, a, a) == kinterval ([2; 2; 2], [5; 5; 5]));
