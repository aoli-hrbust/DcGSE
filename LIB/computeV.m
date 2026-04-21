<<<<<<< HEAD
function V = computeV(U,X)
% [m, n] = size(U);
% [p, q] = size(X);
n = size(U,1);
p = size(X,1);
if n ~= p
    error('矩阵尺寸不匹配。U的列数必须等于X的行数。');
end
V = U' * X;
=======
function V = computeV(U,X)
% [m, n] = size(U);
% [p, q] = size(X);
n = size(U,1);
p = size(X,1);
if n ~= p
    error('矩阵尺寸不匹配。U的列数必须等于X的行数。');
end
V = U' * X;
>>>>>>> e1248d18d20a74f08f96760f80be9bf1df35eab1
end