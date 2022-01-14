## -*- texinfo -*-
## @deftypefn  {} {} kinterval ()
## @deftypefnx {} {} kinterval (@var{a}, @var{b})
## Create a Kaucher interval object
##
## @example
## [a,b]
## @end example
##
## @noindent
## from bounds a and b.
## @end deftypefn

function z = kinterval(varargin)

%  if (nargin > 1)
%    print_usage ();
%  endif

  if (nargin == 0)
    z.inf = 0;
    z.sup = 0;
    z = class (z, "kinterval");
  else
   if nargin == 1
    x = varargin{1}; 
    if (strcmp (class (x), "kinterval"))
      z = x;
    elseif isreal (x)
      z.inf = x;  
      z.sup = x;  
      z = class (z, "kinterval");
    else
        error ("@kinterval: A must be a real numbers");
    endif 
   elseif  nargin == 2 
      x = varargin{1}; 
      y = varargin{2}; 
      if isreal (x) && isreal (y)
        if ~all(size (x) == size (y))
          error("wrong dimensions");
        endif
        z.inf = x;
        z.sup = y;
        z = class (z, "kinterval");
      else
        error ("@kinterval: A and B must be real numbers");
      endif
    else 
      error ("@kinterval: Wrong number of arguments");
    endif    
  endif

endfunction



%classdef kinterval
%  properties
%    kinf = 0;
%    ksup = 0;
%  endproperties
%
%  # TODO
%  # Empty
%  # Entire
%  # Vector/Matrix
%  # Rounding
%  # Construction through strings, etc
%  
%  methods
%    function k = kinterval (a,b)
%      if (nargin == 0)
%          k.kinf = 0;
%          k.ksup = 0;
%      elseif (nargin == 1) && isreal(a)
%          k.kinf = a;
%          k.ksup = a; 
%      elseif (nargin == 2) && isreal(a) && isreal(b)
%            k.kinf = a;
%            k.ksup = b;
%      else
%            error ("@kinterval: A, B must be real numbers");
%      endif
%    endfunction
%
%    function display (k)
%	    printf ("%s = ", inputname (1));
%	    printf ("[%d, %d]\n", k.kinf, k.ksup)
%    endfunction
%
%    function t = in(x, y)
%      if (not (isa (x, "kinterval")))
%        x = kinterval (x);
%      endif
%      if (not (isa (y, "kinterval")))
%        y = kinterval (y);
%      endif
%	    t = x.kinf >= y.kinf && x.ksup <= y.ksup;
%    endfunction
%
%    function t = le(x, y)
%      if (not (isa (x, "kinterval")))
%        x = kinterval (x);
%      endif
%      if (not (isa (y, "kinterval")))
%        y = kinterval (y);
%      endif
%	    t = x.kinf <= y.kinf && x.ksup <= y.ksup;
%    endfunction
%
%    function t = ge(x, y)
%      if (not (isa (x, "kinterval")))
%        x = kinterval (x);
%      endif
%      if (not (isa (y, "kinterval")))
%        y = kinterval (y);
%      endif
%	    t = x.kinf >= y.kinf && x.ksup >= y.ksup;
%    endfunction
%
%    function z = abs(x)
%	    z = max(abs(x.kinf), abs(x.ksup));
%      if (not (isa (x, "kinterval")))
%        x = kinterval (x);
%      endif
%    endfunction
%
%    function z = and(x, y)
%      if (not (isa (x, "kinterval")))
%        x = kinterval (x);
%      endif
%      if (not (isa (y, "kinterval")))
%        y = kinterval (y);
%      endif
%	    z = kinterval (max(x.kinf, y.kinf), min(x.ksup, y.ksup));
%    endfunction
%
%    function z = or(x, y)
%      if (not (isa (x, "kinterval")))
%        x = kinterval (x);
%      endif
%      if (not (isa (y, "kinterval")))
%        y = kinterval (y);
%      endif
%	    z = kinterval (min(x.kinf, y.kinf), max(x.ksup, y.ksup));
%    endfunction
%
%    function z = plus (x, y)
%      if (not (isa (x, "kinterval")))
%        x = kinterval (x);
%      endif
%      if (not (isa (y, "kinterval")))
%        y = kinterval (y);
%      endif
%  
%      l = x.kinf + y.kinf;
%      u = x.ksup + y.ksup;
%
%      z = kinterval (l, u);
%    endfunction
%
%    function z = minus (x, y)
%      if (not (isa (x, "kinterval")))
%        x = kinterval (x);
%      endif
%      if (not (isa (y, "kinterval")))
%        y = kinterval (y);
%      endif
%      
%      l = x.kinf - y.ksup;
%      u = x.ksup - y.kinf;
%
%      #emptyresult = isempty (x) | isempty (y);
%      #l(emptyresult) = inf;
%      #u(emptyresult) = -inf;
%      
%      #l(l == 0) = -0;
%      
%      z = kinterval (l, u);
%    endfunction
%
%    function z = mtimes(x, y)
%      z = x.* y;
%    endfunction
%    
%    function z = times(x, y)
%      if (not (isa (x, "kinterval")))
%        x = kinterval (x);
%      endif
%      if (not (isa (y, "kinterval")))
%        y = kinterval (y);
%      endif
%
%      xP = ispositive(x);
%	    yP = ispositive(y);
%      xN = isnegative(x);
%      yN = isnegative(y);
%      xZ = iszeroin(x);
%      yZ = iszeroin(y);
%      
%      # empty
%      # entire
%
%      if  xP 
%          if yP
%            z = kinterval (x.kinf * y.kinf, x.ksup * y.ksup);
%          elseif yZ 
%            z = kinterval (x.ksup * y.kinf, x.ksup * y.ksup);
%          elseif yN
%            z = kinterval (x.ksup * y.kinf, x.kinf * y.ksup);
%          else
%            z = kinterval (x.kinf * y.kinf, x.kinf * y.ksup);
%          endif
%      elseif xZ
%          if yP
%            z = kinterval (x.kinf * y.ksup, x.ksup * y.ksup);
%          elseif yZ 
%            z = kinterval (min(x.kinf * y.ksup, x.ksup * y.kinf), max(x.kinf * y.kinf, x.ksup * y.ksup));
%          elseif yN
%            z = kinterval (x.ksup * y.kinf, x.kinf * y.kinf);
%          else
%            z = kinterval (0, 0);
%          endif
%      elseif xN
%          if yP
%            z = kinterval (x.kinf * y.ksup, x.ksup * y.kinf);
%          elseif yZ 
%            z = kinterval (x.kinf * y.ksup, x.kinf * y.kinf);
%          elseif yN
%            z = kinterval (x.ksup * y.ksup, x.kinf * y.kinf);
%          else
%            z = kinterval (x.ksup * y.ksup, x.ksup * y.kinf);
%          endif
%      else #if adualZ
%          if yP
%            z = kinterval (x.kinf * y.kinf, x.ksup * y.kinf);
%          elseif yZ 
%            z = kinterval (0, 0);
%          elseif yN
%            z = kinterval (x.ksup * y.ksup, x.kinf * y.ksup);
%          else
%            z = kinterval (max(x.kinf * y.kinf, x.ksup * y.ksup), min(x.kinf * y.ksup, x.ksup * y.kinf));
%          endif
%      endif
%    endfunction
%
%    function z = inv(x)
%      if iszeroin(pro(x))
%        error ('wrong division, pro(denominator) contains zero');
%      else
%        z = kinterval(1/x.kinf, 1/x.ksup);	      
%      endif
%    endfunction
%
%    function z = rdivide(x,y)
%      if (not (isa (x, "kinterval")))
%        x = kinterval (x);
%      endif
%      if (not (isa (y, "kinterval")))
%        y = kinterval (y);
%      endif
%      if iszeroin(pro(y))
%        error ('wrong division, pro(denominator) contains zero');
%      endif
%      inv_y = kinterval(1/y.ksup, 1/y.kinf);
%      z = x .* inv_y;
%    endfunction
%
%    function z = mrdivide(x,y)
%      z = rdivide(x,y);
%    endfunction
%    
%    function z = opp(x)
%	    z = kinterval(-x.kinf, -x.ksup);
%    endfunction
%
%    function z = dual(x)
%	    z = kinterval(x.ksup, x.kinf);
%    endfunction
%
%    function z = pro(x)
%      if x.kinf <= x.ksup
%          z = x;
%      else
%          z = kinterval(x.ksup, x.kinf);
%      endif
%    endfunction
%
%    function t = ispositive(x)
%    	t = x.kinf >= 0 &&  x.ksup >= 0;
%    endfunction
%    
%    function t = iszeroin(x)
%	    t = x.kinf <= 0 &&  0 <= x.ksup;
%    endfunction
%
%    function t = isnegative(x)
%	    t = x.kinf <= 0 &&  x.ksup <= 0;
%    endfunction
%
%    function t = isinzero(x)
%	    t = iszeroin(dual(x));
%    endfunction
%
%  endmethods
%endclassdef