% Jacobian
function j=JAt(xi, eta, r, c)
  global X Y rN cN;
  s = 1+(r-1)*(cN+1)+(c-1);
  j = dNAt(xi, eta)*[X(s) X(s+cN+1) X(s+cN+2) X(s+1);
                     Y(s) Y(s+cN+1) Y(s+cN+2) Y(s+1); ]';
end