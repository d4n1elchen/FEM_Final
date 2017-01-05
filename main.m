clc;
clear;

% plot setting
mesh = false;
dvector = true;

global X Y D T rN cN nN;

rN = 64;              % Number of mesh row
cN = 128;              % Number of mesh coloum
nN = (rN+1)*(cN+1);   % Number of node

a = 0.5;
X = [0 0 2   2];
Y = [1 0 0.5 1];

[X Y] = MESH(X, Y, rN, cN);

nu = 0.3;
E = 3e7;
D = E/(1-nu^2)*[ 1 nu 0;
                 nu 1 0;
                 0  0 (1-nu)/2 ];

T = [ 0; -20 ];

%dNod = [ 1 2 3 4 ];
%dVal = [ 0 0 0 0 ];

mKe = Ke(1, 1);
%F = Fe()

K = Kg(rN, cN);
F = Fg();

d = Inf(nN*2,1);
dNod = [ 1 2 ];
dVal = [ 0 0 ];
d(dNod) = dVal(:);

for i=1:rN
  dNod = [ i*(cN+1)*2+1 i*(cN+1)*2+2 ];
  dVal = [ 0 0 ];
  d(dNod) = dVal(:);
end
d;

res = SOLVE(K, F, d);
[s gp] = STRESS(res);
s = reshape(s,cN*rN*4,3);
gp = reshape(gp,cN*rN*4,2);

r=reshape(res,(cN+1)*2,(rN+1))';
x=r(:,[1:2:end]);
y=r(:,[2:2:end]);
X=reshape(X,(cN+1),(rN+1))';
Y=reshape(Y,(cN+1),(rN+1))';
Z=sqrt(x.^2+y.^2);

row = 2*rN;
col = 2*cN;
figure;

subplot(2,2,1);
[cs, hc] = contourf(X,Y,Z,20);
set(hc,'EdgeColor','none');
colormap(jet);
hold on;
if mesh
  plot(X,Y,'r');
  plot(X',Y','r');
end
if dvector
  quiver(X,Y,x,y);
end
hold off;

subplot(2,2,2);
data = [gp(:,1) gp(:,2) s(:,1)];
data = sortrows(data);
[cs, hc] = contourf(reshape(data(:,1),row,col),reshape(data(:,2),row,col),reshape(data(:,3),row,col),20);
set(hc,'EdgeColor','none');
colormap(jet);
hold on;
if mesh
  plot(X,Y,'r');
  plot(X',Y','r');
end
if dvector
  quiver(X,Y,x,y);
end
hold off;

subplot(2,2,3);
data = [gp(:,1) gp(:,2) s(:,2)];
data = sortrows(data);
[cs, hc] = contourf(reshape(data(:,1),row,col),reshape(data(:,2),row,col),reshape(data(:,3),row,col),20);
set(hc,'EdgeColor','none');
colormap(jet);
hold on;
if mesh
  plot(X,Y,'r');
  plot(X',Y','r');
end
if dvector
  quiver(X,Y,x,y);
end
hold off;

subplot(2,2,4);
data = [gp(:,1) gp(:,2) s(:,3)];
data = sortrows(data);
[cs, hc] = contourf(reshape(data(:,1),row,col),reshape(data(:,2),row,col),reshape(data(:,3),row,col),20);
set(hc,'EdgeColor','none');
colormap(jet);
hold on;
if mesh
  plot(X,Y,'r');
  plot(X',Y','r');
end
if dvector
  quiver(X,Y,x,y);
end
hold off;


% Mesh
function [x y]=MESH(X, Y, rN, cN)
  global nN;
  xx = [linspace(X(1),X(4),cN+1); 
       linspace(X(2),X(3),cN+1);];
  yy = [linspace(Y(1),Y(4),cN+1); 
       linspace(Y(2),Y(3),cN+1);];
  x = [];
  y = [];
  for i=1:cN+1
    x = [x linspace(xx(1, i), xx(2, i), rN+1)'];
    y = [y linspace(yy(1, i), yy(2, i), rN+1)'];
  end
  x = reshape(x',1,nN);
  y = reshape(y',1,nN);
end

% dN/dxi
function b=dNAt(xi, eta)
  b = (1/4)*[eta-1 1-eta 1+eta -eta-1;
             xi-1  -xi-1 1+xi  1-xi;];
end

% Jacobian
function j=JAt(xi, eta, r, c)
  global X Y rN cN;
  s = 1+(r-1)*(cN+1)+(c-1);
  j = dNAt(xi, eta)*[X(s) X(s+cN+1) X(s+cN+2) X(s+1);
                     Y(s) Y(s+cN+1) Y(s+cN+2) Y(s+1); ]';
end

% dN/dx
function r=DNDX(xi, eta, r, c)
  r = inv(JAt(xi, eta, r, c))*dNAt(xi, eta);
end

% N matrix
function n=Ne(xi, eta)
  n = 1/4*[ (1-xi)*(1-eta) 0 (1+xi)*(1-eta) 0 (1+xi)*(1+eta) 0 (1-xi)*(1+eta) 0; 
            0 (1-xi)*(1-eta) 0 (1+xi)*(1-eta) 0 (1+xi)*(1+eta) 0 (1-xi)*(1+eta); ]; 
end

% B matrix
function b=Be(xi, eta, r, c)
  d = DNDX(xi, eta, r, c);
  [m, n] = size(d);
  b = [];
  for i=1:n
    b = [ b [d(1,i) 0;
            0      d(2,i);
            d(2,i) d(1,i);] ];
  end
end

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

% Global K
function k=Kg(rN, cN)
  global nN;
  n = nN;
  k = zeros(n*2);
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

% global force vector
function f=Fg()
  global cN nN;
  f = zeros(nN*2, 1);
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

% Solve
function d=SOLVE(K, F, d)
  mm = [K F];
  i = find(d==0);
  mm(i,:) = [];
  mm(:,i) = [];
  i = find(d==inf);
  d(i) = gaussElim(mm);
end

% Gaus elimination
function x=gaussElim(A)
  [m n] = size(A);
  x = zeros(m,1); 
  for i = 1:m-1
    a = -A(i+1:m,i)/A(i,i);
    A(i+1:m,:) = A(i+1:m,:) + a*A(i,:);
  end;

  x(m) = A(m,n)/A(m,m);
  for i = m-1:-1:1
      x(i) = (A(i,n) - A(i,i+1:m)*x(i+1:m))/A(i,i);
  end
end

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

% TODO: Not yet finish
% xi eta to x y
function [x y]=toXY(xi, eta, X, Y)
  [XI ETA] = meshgrid(-1:2:1);
  x  = interp2(XI, ETA, [X(1) X(2); X(4) X(3);], xi, eta);
  y  = interp2(XI, ETA, [Y(1) Y(2); Y(4) Y(3);], xi, eta);
end

% TODO: Not yet finish
% find point on which element
function [row col]=ELEMOF(x, y)
  
end
