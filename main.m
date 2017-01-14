clc;
clear;

global X Y D T rN cN nN gpu;

% plot setting
mesh = false;
dvector = true;

% gpu support
gpu = false;

rN = 20;              % Number of mesh row
cN = 40;              % Number of mesh coloum
nN = (rN+1)*(cN+1);   % Number of node

a = 0;
X = [0 0 2 2];
Y = [1 0 a 1];

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

if gpu
    gK = gpuArray(K);
    gF = gpuArray(F);
    gd = gpuArray(d);
    res = SOLVE(gK, gF, gd);
else
    res = SOLVE(K, F, d);
end

[ss bp] = STRESS_BOUNDARY(res);
bbp = bp;
sss = ss;
bp = [];
ss = [];
for i=1:cN
  tp = bbp(1,i,:,1);
  ts = sss(1,i,:,:);
  bp = [bp; tp(:)];
  ss = [ss; reshape(ts,21,3)];
end
%{
figure;
plot(bp(:,1), ss(:,1));
title('\sigma_x');
figure;
plot(bp(:,1), ss(:,3));
title('\tau');
%}
figure;
plot(bp(:,1), ss(:,2));
title('\sigma_y');
xlabel('y(m)');
ylabel('\sigma(Pa)');

[s gp] = STRESS_GP(res);
ggp = gp;
ss = s;
gp = [];
s = [];
for i=1:rN
  for j=1:cN
    tp = ggp(i,j,:,:);
    ts = ss(i,j,:,:);
    gp = [gp; reshape(tp,4,2)];
    s = [s; reshape(ts,4,3)];
  end
end

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
title('Displacement');
colorbar;
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
title('\sigma_x');
colorbar;
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
title('\sigma_y');
colorbar;
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
title('\tau_{xy}');
colorbar;
hold off;
