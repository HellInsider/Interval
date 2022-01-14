## -*- texinfo -*-
## @documentencoding UTF-8
## @defop Method {@@kinterval} ne (@var{A}, @var{B})
## @defopx Operator {@@kinterval} {@var{A} ~= @var{B}}
## @defopx Operator {@@kinterval} {@var{A} != @var{B}}
##
## Compare intervals @var{A} and @var{B} for inequality.
##
## True, if @var{A} contains a number which is not a member in @var{B} or vice
## versa.  False, if all numbers from @var{A} are also contained in @var{B} and
## vice versa.
##
## @seealso{@@kinterval/eq, @@kinterval/subset, @@kinterval/interior, @@kinterval/disjoint}
## @end defop

function result = ne (a, b)

  if (nargin ~= 2)
    print_usage ();
    return
  endif

  result = not (eq (a, b));

endfunction

%!assert (ne (kinterval (1, 2), kinterval (1, 3)));
