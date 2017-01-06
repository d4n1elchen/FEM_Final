% global force vector
function f=Fg()
  global cN nN;
  f = zeros(nN*2, 1);
  fprintf('Calc global force vector\n');
  for i=1:cN
    F = Fe(1, i);
    g1 = i;
    gNod = [g1 g1+cN+1 g1+cN+2 g1+1];
    gNod = [2*gNod-1 2*gNod];
    lNod = [1:4];
    lNod = [2*lNod-1 2*lNod];
    f(gNod) = f(gNod) + F(lNod);
  end
end