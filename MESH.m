% Mesh
function [x y]=MESH(X, Y, rN, cN)
  global nN;
  xx = [linspace(X(1),X(4),cN+1); 
       linspace(X(2),X(3),cN+1);];
  yy = [linspace(Y(1),Y(4),cN+1); 
       linspace(Y(2),Y(3),cN+1);];
  x = [];
  y = [];
  fprintf('Meshing\n');
  parfor i=1:cN+1
    x = [x linspace(xx(1, i), xx(2, i), rN+1)'];
    y = [y linspace(yy(1, i), yy(2, i), rN+1)'];
  end
  x = reshape(x',1,nN);
  y = reshape(y',1,nN);
end
