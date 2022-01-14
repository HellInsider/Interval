## -*- texinfo -*-
## @documentencoding UTF-8
## @defop Method {@@kinterval} subsasgn (@var{A}, @var{IDX}, @var{RHS})
## @defopx Operator {@@kinterval} {@var{A}(@var{SUBS}) = @var{RHS}}
##
## Perform the subscripted assignment operation according to the subscript
## specified by @var{IDX}.
##
## The subscript @var{IDX} is expected to be a structure array with fields
## @code{type} and @code{subs}.  Only valid value for @var{type} is
## @code{"()"}.  The @code{subs} field may be either @code{":"} or a cell array
## of index values.
## @seealso{@@kinterval/subsref, @@kinterval/end}
## @end defop

function A = subsasgn (A, S, B)

  if (nargin ~= 3)
    print_usage ();
    return
  endif
  if (not (isa (A, "kinterval")))
    A = kinterval (A);
  endif
  if (not (isa (B, "kinterval")))
    B = kinterval (B);
  endif

  assert (strcmp (S.type, "()"), "only subscripts with parenthesis allowed");


  A.inf = subsasgn (A.inf, S, B.inf);
  A.sup = subsasgn (A.sup, S, B.sup);
# A.inf(A.inf == 0) = -0; # any new elements are stored as [-0, +0] internally

endfunction

%!test
%! A = kinterval (magic (3));
%! A(4, 4) = 42;
%! assert (inf (A), [magic(3),[0;0;0];0,0,0,42]);
%! assert (sup (A), [magic(3),[0;0;0];0,0,0,42]);
