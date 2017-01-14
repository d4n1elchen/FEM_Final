% Stress
function [s bp]=STRESS_BOUNDARY(res)
  global rN cN D X Y gpu;
  s = zeros(1, cN, 4, 3);
  bp = zeros(1, cN, 4, 2);
  if gpu
      s = gpuArray(s);
      bp = gpuArray(bp);
  end
  for i=1:cN
    g1 = i;
    gNod = [g1 g1+cN+1 g1+cN+2 g1+1];
    xx = X(gNod);
    yy = Y(gNod);
    g1 = 1+(i-1)*2;
    gNod = [g1 g1+1 g1+(cN+1)*2 g1+(cN+1)*2+1 g1+(cN+1)*2+2 g1+(cN+1)*2+3 g1+2 g1+3];
    d = res(gNod);
    c = 1;
    xi = -1;
    for k=-1:0.1:1
      eta = k;
      s(1, i, c, :) = D*Be(xi, eta, 1, i)*d;
      [x y] = toXY(xi, eta, xx, yy);
      bp(1, i, c, :) = [x y];
      c = c+1;
    end
  end
end