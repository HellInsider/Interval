## -*- texinfo -*-
## @documentencoding UTF-8
## @defmethod {@@kinterval} cat (@var{DIM}, @var{ARRAY1}, @var{ARRAY2}, @dots{})
##
## Return the concatenation of N-D interval arrays @var{ARRAY1}, @var{ARRAY2},
## … along dimension @var{DIM}.
##
## @example
## @group
## cat (2, kinterval (magic (3)), kinterval (pascal (3)))
##   @result{} 3×6 interval matrix
##
##      [8]   [1]   [6]   [1]   [1]   [1]
##      [3]   [5]   [7]   [1]   [2]   [3]
##      [4]   [9]   [2]   [1]   [3]   [6]
## @end group
## @end example
## @seealso{@@kinterval/horzcat, @@kinterval/vertcat}
## @end defmethod

function result = cat (dim, varargin)

  if (isa (dim, "kinterval"))
    print_usage ();
    return
  endif

  ## Conversion of non-interval parameters to intervals
  interval_idx = cellfun ("isclass", varargin, "kinterval");
  to_convert_idx = not (interval_idx);
  varargin(to_convert_idx) = cellfun (@kinterval, ...
                                      varargin(to_convert_idx), ...
                                      "UniformOutput", false);

  ## Stack intervals along dimension dim
  s = cellfun ("struct", varargin); # struct array
  result = kinterval ();
  result.inf = cat (dim, s.inf);
  result.sup = cat (dim, s.sup);

endfunction

%!assert (size (cat (1, kinterval ([]), kinterval ([]))), [0 0]);
%!assert (cat (1, kinterval (1), kinterval (2)) == kinterval (cat (1, 1, 2)));
%!assert (cat (2, kinterval (1), kinterval (2)) == kinterval (cat (2, 1, 2)));
%!assert (cat (5, kinterval (1), kinterval (2)) == kinterval (cat (5, 1, 2)));
%!assert (cat (1, kinterval (zeros (2, 2, 2)), kinterval (ones (2, 2, 2))) == kinterval (cat (1, zeros (2, 2, 2), ones (2, 2, 2))));
