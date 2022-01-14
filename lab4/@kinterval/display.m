## -*- texinfo -*-
## @deftypefn {} {} display (@var{obj})
## Display the contents of the object @var{obj}.
##
## The Octave interpreter calls the @code{display} function whenever it needs
## to present a class on-screen.  Typically, this would be a statement which
## does not end in a semicolon to suppress output.  For example:
##
## @example
## myobj = myclass (@dots{})
## @end example
##
## User-defined classes should overload the @code{display} method so that
## something useful is printed for a class object.  Otherwise, Octave will
## report only that the object is an instance of its class.
##
## @example
## @group
## myobj = myclass (@dots{})
##   @result{} myobj = <class myclass>
## @end group
## @end example
##
## @seealso{class, subsref, subsasgn}
## @end deftypefn

function display (k)
   if (~strcmp (class (k), "kinterval"))
       disp(k)';
       return;
   endif
     
   kinf =  k.inf;
  [rows,cols] = size(kinf);
  
  printf ("%s = ", inputname (1));

  prec = output_precision();

  if rows*cols == 1
    printf ("\n");
    #printf ("[%d,%d]\n", k.inf, k.sup);
    printf("[%.*g,%.*g]\n", prec, k.inf, prec, k.sup);
  else
    data = zeros(2*cols, rows);
    data(1:2:end,:) = k.inf';
    data(2:2:end,:) = k.sup';

    printf ("\n");
    #for i = 1:rows
    #  printf (" [%f,%f] ", data(:,i));
    #  printf ("\n");
    #endfor
    for i = 1:rows
        for j = 1:cols            
            printf(" [%.*g,%.*g]", prec, k.inf(i,j), prec, k.sup(i,j));
        endfor        
        printf ("\n");
    endfor
    #printf ("]\n");
  endif
endfunction

%!test
%! str = evalc ("x = 1.1; display (x)");
%! assert (str, "x =  1.1000\n");

%!test
%! str = evalc ("display (1.1)");
%! assert (str, " 1.1000\n");

## Test input validation
%!error display ()
%!error display (1,2)
