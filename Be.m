% B matrix
function b=Be(xi, eta, r, c)
  d = DNDX(xi, eta, r, c);
  [m, n] = size(d);
  b = [];
  for i=1:n
    b = [ b [d(1,i) 0;
            0      d(2,i);
            d(2,i) d(1,i);] ];
  end
end