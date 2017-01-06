% xi eta to x y
function [x y]=toXY(xi, eta, X, Y)
  [XI ETA] = meshgrid(-1:2:1);
  x  = interp2(XI, ETA, [X(1) X(2); X(4) X(3);], xi, eta);
  y  = interp2(XI, ETA, [Y(1) Y(2); Y(4) Y(3);], xi, eta);
end