n = 10;

D = diag([50 -50 1 2 3 4 5 6 7 8]);

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
            fprintf(f_eigenVectors, '%.15f\n', eigen(i, j));
            fprintf(f_matrix, '%.15f\n', A(i, j));
        else
            fprintf(f_eigenVectors, '%.15f ', eigen(i, j));
            fprintf(f_matrix, '%.15f ', A(i, j));
        end
    end
    if (i == n)
        fprintf(f_eigenValues, '%.15f\n', S(i, i));
    else
        fprintf(f_eigenValues, '%.15f ', S(i, i));
    end
end
fclose(f_matrix);
fclose(f_eigenVectors);
fclose(f_eigenValues);