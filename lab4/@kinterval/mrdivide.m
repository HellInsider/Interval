function z = mrdivide(x, y)
 
  if (not (isa (x, "kinterval")))
    x = kinterval (x);
  endif
  if (not (isa (y, "kinterval")))
    y = kinterval (y);
  endif

  if (isscalar (x) || isscalar (y))
      z = x ./ y;
      return
  endif

  [xrows, xcols] = size(x.inf);
  [yrows, ycols] = size(y.inf);
  if  ndims(x.inf) > 2 || ndims(y.inf) > 2 || xcols ~= yrows 
    error("wrong dimensions");
  endif

  error('matrix right division is not implemented yet');

#      zinf = zeros(xrows,ycols);
#      zsup = zeros(xrows,ycols);
#      
#      for i = 1:xrows
#          for j = 1: ycols
#              p = x(i,:) .* y(:,j)';
#              zinf(i,j) = sum(p.inf);
#              zsup(i,j) = sum(p.sup);              
#          endfor
#      endfor
#                
#      z = kinterval(zinf, zsup);

endfunction
