n = 3;

D = diag([7 -7 1]);

[Q, R] = qr(rand(n));
A = Q*D*Q';
[eigen, S] = eig(A);
eigen = eigen';
f_matrix = fopen('matrix.txt', 'w');
f_eigenVectors = fopen('eigenVectorsKnown.txt', 'w');
f_eigenValues = fopen('eigenValuesKnown.txt', 'w');
for i=1:n
    for j=1:n
        if (j == n)
            fprintf(f_eigenVectors, '%.4f\n', eigen(i, j));
            fprintf(f_matrix, '%.4f\n', A(i, j));
        else
            fprintf(f_eigenVectors, '%.4f ', eigen(i, j));
            fprintf(f_matrix, '%.4f ', A(i, j));
        end
    end
    if (i == n)
        fprintf(f_eigenValues, '%.4f\n', S(i, i));
    else
        fprintf(f_eigenValues, '%.4f ', S(i, i));
    end
end
fclose(f_matrix);
fclose(f_eigenVectors);
fclose(f_eigenValues);