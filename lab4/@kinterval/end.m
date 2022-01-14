## -*- texinfo -*-
## @documentencoding UTF-8
## @defmethod {@@kinterval} end (@var{A}, @var{k}, @var{n})
##
## The magic index @code{end} refers to the last valid entry in an indexing
## operation.
##
## For @var{n} = 1 the output of the indexing operation is a column vector and a
## single index is used to address all entries in column-first order.  For
## @var{n} > 1 the @var{k}'th dimension is addressed separately.
##
## @example
## @group
## A = kinterval ([1 2 3; 4 5 6]);
## A(1, end)
##   @result{} ans =  [3]
## A(end, 1)
##   @result{} ans =  [4]
## A(end)
##   @result{} ans =  [6]
## @end group
## @end example
## @seealso{@@kinterval/size, @@kinterval/length, @@kinterval/numel, @@kinterval/rows, @@kinterval/columns}
## @end defmethod

function ret = end (a, k, n)

  if (n == k)
    ret = prod (size (a.inf)(n:ndims (a.inf)));
  else
    ret = size (a.inf, k);
  endif

endfunction

%!assert (kinterval (magic (3))(end) == 2);
%!assert (kinterval (magic (3))(end, 2) == 9);
%!assert (kinterval (magic (3))(2, end) == 7);
%!assert (kinterval ([1 2; 3 4; 5 6])(end:-1:1, :) == [5 6; 3 4; 1 2]);
%!assert (reshape (kinterval (1:24), 2, 3, 4)(end, end) == 24)
