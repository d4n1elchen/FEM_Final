% Solve
function d=SOLVE(K, F, d)
  mm = [K F];
  i = find(d==0);
  mm(i,:) = [];
  mm(:,i) = [];
  i = find(d==inf);
  fprintf('Solve the equation\n');
  d(i) = gaussElim(mm);
end