function z = dist(x,y)
    if (not (isa (x, "kinterval")))
        x = kinterval (x);
    endif
    if (not (isa (y, "kinterval")))
        y = kinterval (y);
    endif
    z = abs(innerminus(x,y));
endfunction

%!# from the documentation string
%!assert (dist (3.6, kinterval (3, 5)) == 1.4);
%!assert (dist (kinterval (3, 5), kinterval (5,3)) == 2);
