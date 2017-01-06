% Stress
function [s gp]=STRESS(res)
  global rN cN nN D X Y;
  g = [ -1/sqrt(3) 1/sqrt(3) ];
  s = zeros(rN, cN, 4, 3);
  gp = zeros(rN, cN, 4, 2);
  for i=1:rN
    for j=1:cN
      g1 = (i-1)*(cN+1)+j;
      gNod = [g1 g1+cN+1 g1+cN+2 g1+1];
      xx = X(gNod);
      yy = Y(gNod);
      g1 = 1+(i-1)*(cN+1)*2+(j-1)*2;
      gNod = [g1 g1+1 g1+(cN+1)*2 g1+(cN+1)*2+1 g1+(cN+1)*2+2 g1+(cN+1)*2+3 g1+2 g1+3];
      d = res(gNod);
      c = 1;
      for k=1:2
        for l=1:2
          xi  = g(k);
          eta = g(l);
          s(i, j, c, :) = D*Be(xi, eta, i, j)*d;
          [x y] = toXY(xi, eta, xx, yy);
          gp(i, j, c, :) = [x y];
          c = c + 1;
        end
      end
    end
  end
end