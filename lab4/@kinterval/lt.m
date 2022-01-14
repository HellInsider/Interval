## -*- texinfo -*-
## @documentencoding UTF-8
## @defop Method {@@kinterval} lt (@var{A}, @var{B})
## @defopx Operator {@@kinterval} {@var{A} < @var{B}}
##
## Compare intervals @var{A} and @var{B} for strictly less.
##
## True, if all numbers from @var{A} are strict less than at least one number
## in @var{B}.  False, if @var{A} contains a number which is greater than all
## numbers in @var{B} or is equal to the greatest number of @var{B}.
##
## Evaluated on interval arrays, this functions is applied element-wise.
##
## @seealso{@@kinterval/eq, @@kinterval/le, @@kinterval/gt, @@kinterval/subset, @@kinterval/interior, @@kinterval/disjoint}
## @end defop

function result = lt (a, b)

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

  result = ((a.inf < b.inf | (a.inf == -inf & b.inf == -inf)) & ...
            (a.sup < b.sup | (a.sup == inf & b.sup == inf)));

#  result (isempty (a) & isempty (b)) = true ();

endfunction

%!assert (not (lt (kinterval (1, 3), kinterval (3))));
%!assert (lt (kinterval (1, 3), kinterval (3.1)));

%!shared testdata
%! # Load compiled test data (from src/test/*.itl)
%! testdata = load (file_in_loadpath ("test/itl.mat"));

%!test
%! # Scalar evaluation
%! testcases = testdata.NoSignal.kinterval.strictLess;
%! for testcase = [testcases]'
%!   assert (isequaln (...
%!     lt (testcase.in{1}, testcase.in{2}), ...
%!     testcase.out));
%! endfor

%!test
%! # Vector evaluation
%! testcases = testdata.NoSignal.kinterval.strictLess;
%! in1 = vertcat (vertcat (testcases.in){:, 1});
%! in2 = vertcat (vertcat (testcases.in){:, 2});
%! out = vertcat (testcases.out);
%! assert (isequaln (lt (in1, in2), out));

%!test
%! # N-dimensional array evaluation
%! testcases = testdata.NoSignal.kinterval.strictLess;
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
%! assert (isequaln (lt (in1, in2), out));
