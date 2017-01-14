% Gaus elimination
function x=gaussElim(A)
  global gpu;
  [m n] = size(A);
  x = zeros(m,1);
  if gpu
    x = gpuArray(x);
  end
  fprintf('Gaus elimination\n');
  textprogressbar('Gaus elimination: ');
  tic;
  for i = 1:m-1
    a = -A(i+1:m,i)/A(i,i);
    A(i+1:m,:) = A(i+1:m,:) + a*A(i,:);
    textprogressbar(i/(m-1)*100);
  end;
  textprogressbar('done');
  toc

  fprintf('Find unknown\n');
  textprogressbar('Find unknown: ');
  tic;
  x(m) = A(m,n)/A(m,m);
  for i = m-1:-1:1
    x(i) = (A(i,n) - A(i,i+1:m)*x(i+1:m))/A(i,i);
    textprogressbar((m-i)/(m-1)*100);
  end
  textprogressbar('done');
  toc
end