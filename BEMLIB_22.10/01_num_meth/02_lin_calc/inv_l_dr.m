%=======================================
% Inverse of a lower triangular L matrix
%=======================================

%----------------
% read the matrix
%----------------
                                                                                
file1 = fopen('matrix_l.dat');
                                                                                
N = fscanf(file1,'%f',[1,1]);
                                                                                
for i=1:N
 for j=1:i
  a(i,j) = fscanf(file1,'%f',[1,1]);
 end
end
                                                                                
fclose(file1);

%-----
% zero the upper triangular
%----
                                                                                
for i=1:N
 for j=i+1:N
  a(i,j) = 0;
 end
end
                                                                                
a

%--------------------
% compute the inverse
%--------------------
                                                                                
b = inv_l(N,a)
                                                                                
%-------
% verify ab=I
%-------
                                                                                
for i=1:N
  for j=1:N
    ab(i,j) = 0.0;
      for  k=1:N
        ab(i,j) = ab(i,j) + a(i,k)*b(k,j);
      end
  end
end

ab
