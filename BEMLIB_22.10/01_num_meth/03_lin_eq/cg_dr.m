clear all

%---
% driver for conjugate gradients
%---

file1 = fopen('mat_s_vec.dat');
N   = fscanf(file1,'%f',[1,1]);
A   = fscanf(file1,'%f',[N,N]);
rhs = fscanf(file1,'%f',[1,N]);
fclose(file1);

%---
% solve
%---

[solution] = cg (N, A, rhs);

%---
% print and confirm
%---

solution'
disp ('Ax-b:'); A*solution'-rhs'
