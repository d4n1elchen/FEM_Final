% Global K
function k=Kg(rN, cN)
  global nN;
  n = nN;
  k = zeros(n*2);
  fprintf('Calc global K\n');
  for i=1:rN
    for j=1:cN
      ke = Ke(i, j);
      g1 = 1+(i-1)*(cN+1)+(j-1);
      gNod = [g1 g1+cN+1 g1+cN+2 g1+1];
      gNod = [2*gNod-1 2*gNod];
      lNod = [1:4];
      lNod = [2*lNod-1 2*lNod];
      k(gNod, gNod) = k(gNod, gNod) + ke(lNod, lNod);
    end
  end
end