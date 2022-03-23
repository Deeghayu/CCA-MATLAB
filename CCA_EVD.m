function [U,V,rho] = CCA_EVD(R11,R12,R22)

R21 = transpose(R12);

A = R11^(-1/2)*R12*R22^(-1/2);
B = A'*A;
C = A*A';

% E = R11^(-1/2)*R12*R22^(-1)*R21*R11^(-1/2);
% C = E

[V,lamda_V] = eig(C);
eigen_valsV = diag(lamda_V);
[eigen_valsV,ind] = sort(eigen_valsV,'descend');
V = V(:,ind);
lamda_V_new = diag(eigen_valsV);

[U,lamda_U] = eig(B);
eigen_valsU = diag(lamda_U);
[eigen_valsU,ind] = sort(eigen_valsU,'descend');
U = U(:,ind);
lamda_U_new = diag(eigen_valsU);

rho = sqrt(lamda_V_new);
% U = R11^(-1/2)*U;
% V = R22^(-1/2)*V;
end

