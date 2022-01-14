function z = mtimes(x, y)

  if (not (isa (x, "kinterval")))
    x = kinterval (x);
  endif
  if (not (isa (y, "kinterval")))
    y = kinterval (y);
  endif

  if (isscalar (x) || isscalar (y))
      z = times (x, y);
      return
  endif

  [xrows, xcols] = size(x.inf);
  [yrows, ycols] = size(y.inf);
  if  ndims(x.inf) > 2 || ndims(y.inf) > 2 || xcols ~= yrows
    error("wrong dimensions");
  endif

  zinf = zeros(xrows,ycols);
  zsup = zinf;

  for i = 1:xrows
      for j = 1: ycols
          #p = x(i,:) .* y(:,j)';

          #r = x(i,:);
          ridx.type = "()"; 
          ridx.subs = {i,":"};
          r = subsref(x,ridx);

          #c = y(:,j)';
          cidx.type = "()"; 
          cidx.subs = {":",j};
          c = subsref(y,cidx);

          p = r .* c';

          zinf(i,j) = sum(p.inf);
          zsup(i,j) = sum(p.sup);
      endfor
  endfor

  z = kinterval(zinf, zsup);
endfunction

%    % Assumes dimensions are compatible
%    function z = dotproduct(x,y)
%        z = sumvec(x.*y);        
%    endfunction
%    
%    % Assumes dimensions are compatible
%    function z = sumvec(x)
%        z = kinterval(sum(x.inf), sum(x.sup));        
%    endfunction
