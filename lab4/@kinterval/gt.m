## -*- texinfo -*-
## @documentencoding UTF-8
## @defop Method {@@kinterval} gt (@var{A}, @var{B})
## @defopx Operator {@@kinterval} {@var{A} > @var{B}}
##
## Compare intervals @var{A} and @var{B} for strict greater.
##
## Evaluated on interval arrays, this functions is applied element-wise.
##
## @seealso{@@kinterval/eq, @@kinterval/lt, @@kinterval/ge}
## @end defop

function result = gt (a, b)

  if (nargin ~= 2)
    print_usage ();
    return
  endif

  result = lt (b, a);

endfunction

%!assert (not (gt (kinterval (2, 3), kinterval (1, 3))));
%!assert (gt (kinterval (2, 3.1), kinterval (1, 3)));
