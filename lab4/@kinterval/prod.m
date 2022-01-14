## -*- texinfo -*-
## @documentencoding UTF-8
## @defmethod {@@kinterval} prod (@var{X})
## @defmethodx {@@kinterval} prod (@var{X}, @var{DIM})
##
## Product of elements along dimension @var{DIM}.  If @var{DIM} is omitted, it
## defaults to the first non-singleton dimension.
##
## Accuracy: The result is a valid enclosure.
##
## @example
## @group
## prod (kinterval (1 : 4))
##   @result{} ans = [24]
## @end group
## @end example
## @seealso{@@kinterval/sum}
## @end defmethod

function [result] = prod (x, dim)

  if (nargin > 2)
    print_usage ();
    return
  endif

  if (nargin < 2)
    ## Try to find non-singleton dimension
    dim = find (size (x.inf) ~= 1, 1);
    if (isempty (dim))
      dim = 1;
    endif
  endif

  if (dim > ndims (x))
    ## Nothing to do
    result = x;
    return
  endif

  if (ndims (x) == 2 && size (x) == [0 0])
    ## Inconsistency: prod ([]) = 1
    x = kinterval(ones (0, 1));
  endif

  resultsize = size (x);
  resultsize (dim) = 1;
  result = kinterval (ones (resultsize));

  idx.type = "()";
  idx.subs = repmat ({":"}, ndims (x));
  for i = 1:size (x, dim)
    idx.subs{dim} = i;
    result = result .* subsref (x, idx);
  endfor

endfunction

%!# from the documentation string
%!assert (prod (kinterval (1 : 4)) == 24);

%!assert (prod (kinterval ([])) == 1);

%!assert (prod (kinterval (magic (3))) == [96, 45, 84]);
%!assert (prod (kinterval (magic (3)), 2) == [48; 105; 72]);
%!assert (prod (kinterval (magic (3)), 3) == magic (3));

%!assert (prod (prod (reshape (kinterval (1:24), 1, 2, 3, 4))) == reshape ([720, 665280, 13366080, 96909120], 1, 1, 1, 4))
