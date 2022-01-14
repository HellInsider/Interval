## -*- texinfo -*-
## @documentencoding UTF-8
## @defop Method {@@infsup} eq (@var{A}, @var{B})
## @defopx Operator {@@infsup} {@var{A} == @var{B}}
##
## Compare intervals @var{A} and @var{B} for equality.
##
## True, if all numbers from @var{A} are also contained in @var{B} and vice
## versa.  False, if @var{A} contains a number which is not a member in @var{B}
## or vice versa.
##
## @example
## @group
## x = infsup (2, 3);
## y = infsup (1, 2);
## x == y
##   @result{} ans =  0
## @end group
## @end example
## @seealso{@@infsup/subset, @@infsup/interior, @@infsup/disjoint}
## @end defop

function result = eq (a, b)

  if (nargin ~= 2)
    print_usage ();
    return
  endif
  if (not (isa (a, "kinterval")))
    a = kinterval (a);
  endif
  if (not (isa (b, "kinterval")))
    b = kinterval (b);
  endif

  result = (a.inf == b.inf & a.sup == b.sup);

endfunction

%!test "Empty interval";
%! assert (eq (kinterval (), kinterval ()) == true);
%! assert (eq (kinterval (), kinterval (1)) == false);
%! assert (eq (kinterval (0), kinterval ()) == false);
%! assert (eq (kinterval (-inf, inf), kinterval ()) == false);
%!test "Singleton intervals";
%! assert (eq (kinterval (0), kinterval (1)) == false);
%! assert (eq (kinterval (0), kinterval (0)) == true);
%!test "Bounded intervals";
%! assert (eq (kinterval (1, 2), kinterval (3, 4)) == false);
%! assert (eq (kinterval (1, 2), kinterval (2, 3)) == false);
%! assert (eq (kinterval (1, 2), kinterval (1.5, 2.5)) == false);
%! assert (eq (kinterval (1, 2), kinterval (1, 2)) == true);
%!test "Unbounded intervals";
%! assert (eq (kinterval (0, inf), kinterval (-inf, 0)) == false);
%! assert (eq (kinterval (0, inf), kinterval (0, inf)) == true);
%! assert (eq (kinterval (-inf, 0), kinterval (-inf, 0)) == true);
%! assert (eq (kinterval (-inf, inf), kinterval (42)) == false);
%! assert (eq (kinterval (-inf, 0), kinterval (-inf, inf)) == false);
%! assert (eq (kinterval (-inf, inf), kinterval (-inf, inf)) == true);

%!shared testdata
%! # Load compiled test data (from src/test/*.itl)
%! testdata = load (file_in_loadpath ("test/itl.mat"));

%!test
%! # Scalar evaluation
%! testcases = testdata.NoSignal.kinterval.equal;
%! for testcase = [testcases]'
%!   assert (isequaln (...
%!     eq (testcase.in{1}, testcase.in{2}), ...
%!     testcase.out));
%! endfor

%!test
%! # Vector evaluation
%! testcases = testdata.NoSignal.kinterval.equal;
%! in1 = vertcat (vertcat (testcases.in){:, 1});
%! in2 = vertcat (vertcat (testcases.in){:, 2});
%! out = vertcat (testcases.out);
%! assert (isequaln (eq (in1, in2), out));

%!test
%! # N-dimensional array evaluations
%! testcases = testdata.NoSignal.kinterval.equal;
%! in1 = vertcat (vertcat (testcases.in){:, 1});
%! in2 = vertcat (vertcat (testcases.in){:, 2});
%! out = vertcat (testcases.out);
%! # Reshape data
%! i = -1;
%! do
%!   i = i + 1;
%!   testsize = factor (numel (in1) + i);
%! until (numel (testsize) > 2)
%! in1 = reshape ([in1; in1(1:i)], testsize);
%! in2 = reshape ([in2; in2(1:i)], testsize);
%! out = reshape ([out; out(1:i)], testsize);
%! assert (isequaln (eq (in1, in2), out));
