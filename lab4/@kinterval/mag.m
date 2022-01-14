## -*- texinfo -*-
## @documentencoding UTF-8
## @defmethod {@@kinterval} mag (@var{X})
##
## Get the magnitude of numbers in interval @var{X}, that is the maximum of
## absolute values for each element.
##
## If @var{X} is empty, @code{mag (@var{X})} is NaN.
##
## Accuracy: The result is exact.
##
## @example
## @group
## mag (kinterval (-3, 2))
##   @result{} ans = 3
## @end group
## @end example
## @seealso{@@kinterval/mig, @@kinterval/smig, @@kinterval/inf, @@kinterval/sup}
## @end defmethod

function result = mag (x)

  if (nargin ~= 1)
    print_usage ();
    return
  endif

  result = max (abs (x.inf), abs (x.sup));
  #result(x.inf > x.sup) = nan ();

endfunction

%!assert (mag (kinterval (-3, 4)), 4);
%!assert (mag (kinterval (-4, 3)), 4);

%!shared testdata
%! # Load compiled test data (from src/test/*.itl)
%! testdata = load (file_in_loadpath ("test/itl.mat"));

%!test
%! # Scalar evaluation
%! testcases = testdata.NoSignal.kinterval.mag;
%! for testcase = [testcases]'
%!   assert (isequaln (...
%!     mag (testcase.in{1}), ...
%!     testcase.out));
%! endfor

%!test
%! # Vector evaluation
%! testcases = testdata.NoSignal.kinterval.mag;
%! in1 = vertcat (vertcat (testcases.in){:, 1});
%! out = vertcat (testcases.out);
%! assert (isequaln (mag (in1), out));

%! # N-dimensional array evaluation
%! testcases = testdata.NoSignal.kinterval.mag;
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
%! assert (isequaln (mag (in1), out));

%!test
%! # Decorated N-dimensional array evaluation
%! testcases = testdata.NoSignal.kinterval.mag;
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
%! assert (isequaln (mag (in1), out));
