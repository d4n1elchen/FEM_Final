% Gaus elimination
function x=gaussElim(A)
  global gpu;
  [m n] = size(A);
  x = zeros(m,1);
  if gpu
    x = gpuArray(x);
  end
  fprintf('Gaus elimination\n');
  for i = 1:m-1
    a = -A(i+1:m,i)/A(i,i);
    A(i+1:m,:) = A(i+1:m,:) + a*A(i,:);
  end;

  fprintf('Find unknown\n');
  x(m) = A(m,n)/A(m,m);
  for i = m-1:-1:1
    x(i) = (A(i,n) - A(i,i+1:m)*x(i+1:m))/A(i,i);
  end
end