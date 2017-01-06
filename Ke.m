% local K
function k=Ke(r, c)
  global D;
  g = [ -1/sqrt(3) 1/sqrt(3) ];
  k = zeros(8, 8);
  for i=1:2
    for j=1:2
      xi  = g(i);
      eta = g(j);
      k = k + Be(xi, eta, r, c)'*D*Be(xi, eta, r, c)*det(JAt(xi, eta, r, c));
    end
  end
end