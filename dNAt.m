% dN/dxi
function b=dNAt(xi, eta)
  b = (1/4)*[eta-1 1-eta 1+eta -eta-1;
             xi-1  -xi-1 1+xi  1-xi;];
end