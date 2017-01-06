% dN/dx
function r=DNDX(xi, eta, r, c)
  r = inv(JAt(xi, eta, r, c))*dNAt(xi, eta);
end