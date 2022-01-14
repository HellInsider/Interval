## -*- texinfo -*-
## @documentencoding UTF-8
## @deftypemethod {@@kinterval} {@var{R} =} rad (@var{X})
## @deftypemethodx {@@kinterval} {[@var{M}, @var{R}] =} rad (@var{X})
##
## Get the radius (and midpoint) of interval @var{X}.
##
## If @var{X} is empty, @code{rad (@var{X})} is NaN.
## If @var{X} is unbounded in one or both directions, @code{rad (@var{X})} is
## positive infinity.
##
## With two output parameters the midpoint and radius are computed.
##
## Accuracy: The radius @var{R} will make a tight enclosure of the interval
## together with midpoint @var{M}.
##
## @example
## @group
## rad (kinterval (2.5, 3.5))
##   @result{} ans = 0.50000
## @end group
## @end example
## @seealso{@@kinterval/inf, @@kinterval/sup, @@kinterval/mid, @@kinterval/wid}
## @end deftypemethod

function [out1, out2] = rad (x)

  if (nargin ~= 1)
    print_usage ();
    return
  endif

  m = mid (x);
  ## The midpoint is rounded to nearest and the radius must cover both boundaries
  #r1 = mpfr_function_d ('minus', +inf, m, x.inf);
  #r2 = mpfr_function_d ('minus', +inf, x.sup, m);
  r1 = abs(m - x.inf);
  r2 = abs(x.sup - m);
  r = max (r1, r2);

  #r (isempty (x)) = nan ();

  if (nargout >= 2)
    out1 = m;
    out2 = r;
  else
    out1 = r;
  endif

endfunction

%!test;
%! [m, r] = rad (kinterval (2.5, 3.5));
%! assert (m, 3);
%! assert (r, .5);

%!# from the documentation string
%!assert (rad (kinterval (2.5, 3.5)), .5);

%!shared testdata
%! # Load compiled test data (from src/test/*.itl)
%! testdata = load (file_in_loadpath ("test/itl.mat"));

%!test
%! # Scalar evaluation
%! testcases = testdata.NoSignal.kinterval.rad;
%! for testcase = [testcases]'
%!   assert (isequaln (...
%!     rad (testcase.in{1}), ...
%!     testcase.out));
%! endfor

%!test
%! # Vector evaluation
%! testcases = testdata.NoSignal.kinterval.rad;
%! in1 = vertcat (vertcat (testcases.in){:, 1});
%! out = vertcat (testcases.out);
%! assert (isequaln (rad (in1), out));

%!test
%! # N-dimensional array evaluation
%! testcases = testdata.NoSignal.kinterval.rad;
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
%! assert (isequaln (rad (in1), out));

%!test
%! # Decorated N-dimensional array evaluation
%! testcases = testdata.NoSignal.kinterval.rad;
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
%! assert (isequaln (rad (in1), out));
