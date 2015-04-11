
n=20;
X = rand(n,1);

xcol = X;
[Xd, ord] = sort(X, 'ascend');
Xdrow = reshape(Xd(1:(n-1)), 1, n-1);
pieceI = repmat(xcol, 1, n-1);
pieceJ = repmat(Xdrow, n, 1);
Delta_matrix = max(pieceI - pieceJ, 0);
Delta_matrix = Delta_matrix(ord,:);
Delta2_matrix = Delta_matrix - (1/n)*repmat(sum(Delta_matrix, 1), n, 1);


M = Delta_matrix';
M = M(2:(n-1), 1:n);

tmp = M'*inv(M*M');


Dsq = zeros(n-2, n);
gaps = Xd(2:n) - Xd(1:(n-1));
%gaps = ones(n-1,1);
%gaps(2) = 0.5;

for ir = 1:(n-2) 
    Dsq(ir, :) = [zeros(1,ir-1), -1/gaps(ir), ...
                    1/gaps(ir) + 1/gaps(ir+1), ...
                    -1/gaps(ir+1), ...
                    zeros(1, n-(ir-1)-3)];
    %Dsq(ir,:) = Dsq(ir,:)/(1/gaps(ir)+1/gaps(ir+1));
end

%Dsq is (n-2)--by--n

tmp = Dsq';
tmp = tmp(3:(n-2), 2:(n-3));

tmp = Delta_matrix(1:n, 2:(n-1));

Dtrunc = [1/gaps(2) + 1/gaps(3), -1/gaps(3), 0;
          -1/gaps(3), 1/gaps(3) + 1/gaps(4), -1/gaps(4);
          0, -1/gaps(4), 1/gaps(4) + 1/gaps(5)];
      

Dtr = Dsq';
Dtr = Dtr(2:(n-1),:);

A = zeros(n-2,n-2);
for i=1:(n-2)
    for j=i:(n-2)
        A(i,j) = sum(gaps(1:i))*sum(gaps((j+1):(n-1)));
    end
end

A = A+A';
for i=1:(n-2)
    A(i,i) = A(i,i)/2;
end
      
      