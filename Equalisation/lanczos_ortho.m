function [] = lanczos_ortho(A, m, debug)

[n,k] = size(A);
V = zeros(k,m+1);
if (nargin == 2)
   V(:,2) = rand(k,1);
else
   V(:,2)=  0.5*ones(k,1);
end
V(:,2)=V(:,2)/norm(V(:,2),2);
beta(2)=0;

for j=2:m+2

    w = A*V(:,j);
    w = A'*w;
    w = w - beta(j)*V(:,j-1);
    alpha(j) = w'*V(:,j);
    w = w - alpha(j)*V(:,j);

    %orthogonalize
    for k=2:j-1
      tmpalpha = w'*V(:,k);
      w = w -tmpalpha*V(:,k);
    end

    beta(j+1) = norm(w,2);
    V(:,j+1) = w/beta(j+1);
end
% prepare the tridiagonal matrix
T = diag(alpha(2:end)) + diag(beta(3:end-1),1) + diag(beta(3:end-1),-1);
V = V(:,2:end-1);
disp(['approximation quality is: ', num2str(norm(V*T*V'-A'*A))]);
disp(['approximating eigenvalues are: ', num2str(sqrt(eig(full(T))'))]);
eigens = eig(full(T));
[u,d]=eig(full(T));
V
disp(['ritz eigenvectors are: ']);
s=V*u
disp(['and u is: ']);
u
disp(['d is: ']);
d
for i=m+1:-1:1
  disp(['residual for ', num2str(i), ' is: ', num2str(norm(A'*A*s(:,i) - eigens(i)*s(:,i)))]);
end