## -*- texinfo -*-
## @documentencoding UTF-8
## @defmethod {@@kinterval} wid (@var{X})
##
## Get the width of interval @var{X}.
##
## If @var{X} is empty, @code{wid (@var{X})} is NaN.
## If @var{X} is unbounded in one or both directions, @code{wid (@var{X})} is
## positive infinity.
##
## Accuracy: The result is a tight enclosure of the interval's actual width.
##
## @example
## @group
## wid (kinterval (2.5, 3.5))
##   @result{} ans = 1
## @end group
## @end example
## @seealso{@@kinterval/inf, @@kinterval/sup, @@kinterval/rad}
## @end defmethod

function result = wid (x)

  if (nargin ~= 1)
    print_usage ();
    return
  endif

  #result = mpfr_function_d ('minus', +inf, x.sup, x.inf);
  #result(x.inf > x.sup) = nan ();

  result = abs(x.sup - x.inf);

endfunction

%!# from the documentation string
%!assert (wid (kinterval (2.5, 3.5)), 1);

%!shared testdata
%! # Load compiled test data (from src/test/*.itl)
%! testdata = load (file_in_loadpath ("test/itl.mat"));

%!test
%! # Scalar evaluation
%! testcases = testdata.NoSignal.kinterval.wid;
%! for testcase = [testcases]'
%!   assert (isequaln (...
%!     wid (testcase.in{1}), ...
%!     testcase.out));
%! endfor

%!test
%! # Vector evaluation
%! testcases = testdata.NoSignal.kinterval.wid;
%! in1 = vertcat (vertcat (testcases.in){:, 1});
%! out = vertcat (testcases.out);
%! assert (isequaln (wid (in1), out));

%!test
%! # N-dimensional array evaluation
%! testcases = testdata.NoSignal.kinterval.wid;
%! in1 = vertcat (vertcat (testcases.in){:, 1});
%! out = vertcat (testcases.out);
%! # Reshape data
%! i = -1;
%! do
%!   i = i + 1;
%!   testsize = factor (numel (in1) + i);
%! until (numel (testsize) > 2)
%! in1 = reshape ([in1; in1(1:i)], testsize);
%! out = reshape ([out; out(1:i)], testsize);
%! assert (isequaln (wid (in1), out));

%!test
%! # Decorated N-dimensional array evaluation
%! testcases = testdata.NoSignal.kinterval.wid;
%! in1 = vertcat (vertcat (testcases.in){:, 1});
%! out = vertcat (testcases.out);
%! # Reshape data
%! i = -1;
%! do
%!   i = i + 1;
%!   testsize = factor (numel (in1) + i);
%! until (numel (testsize) > 2)
%! in1 = reshape ([in1; in1(1:i)], testsize);
%! out = reshape ([out; out(1:i)], testsize);
%! assert (isequaln (wid (in1), out));
