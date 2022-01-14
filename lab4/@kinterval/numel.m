## -*- texinfo -*-
## @documentencoding UTF-8
## @defmethod {@@kinterval} numel (@var{A})
##
## Return the number of elements in the interval object @var{A}.
## @seealso{@@kinterval/length, @@kinterval/size, @@kinterval/rows, @@kinterval/columns, @@kinterval/end}
## @end defmethod

function result = numel (a)

  if (nargin ~= 1)
    print_usage ();
    return
  endif

  result = numel (a.inf);

endfunction

%!assert (numel (kinterval ([])), 0);
%!assert (numel (kinterval (0)), 1);
%!assert (numel (kinterval (zeros (3, 1))), 3);
%!assert (numel (kinterval (zeros (1, 4))), 4);
%!assert (numel (kinterval (zeros (3, 4))), 12);
