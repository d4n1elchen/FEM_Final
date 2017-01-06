% local force vector
function f=Fe(r, c)
  global T;
  xi = -1;
  g = [ -1/sqrt(3) 1/sqrt(3) ];
  w = [ 1 1 ];
  f = zeros(8, 1);
  [m n] = size(g);
  for i=1:n
    eta = g(i);
    J = JAt(xi, eta, r, c);
    f = f + w(i)*Ne(xi, eta)'*T*det(J(2, 1));
  end
end