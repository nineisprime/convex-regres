n_ls = 100:50:500;

s = 3;
p = 40;
kappa_ls = zeros(length(n_ls),1);

num_trial = 10000;

for ii=1:length(n_ls)
    cur_n = n_ls(ii);
    
    X = randn(p,cur_n);
    
    %generate Bstar
    Bstar_candidates = randn(3*cur_n, s);
    Bstar_candidates = [Bstar_candidates,zeros(3*cur_n,p-s)];
    total_prods = Bstar_candidates*X;
    [vals,ixs] = max(total_prods,[],1);
    Bstar = Bstar_candidates(ixs,:);
    %Bstar = zeros(cur_n, p);
    
    kappa = 10e10;
    for it = 1:num_trial        
        B_candidates = randn(3*cur_n,s);
        B_candidates = [B_candidates, zeros(3*cur_n,p-s)];
        B_candidates(:,(s+1):p) = randn(3*cur_n, p-s)*0.4;
        total_prods = B_candidates*X;
        [vals,ixs] = max(total_prods,[],1);
        B = B_candidates(ixs,:);
        
        LHS = (1/cur_n)*sum(diag((B-Bstar)*X).^2);
        %RHS = norm(B-Bstar,'fro')^2;
        RHS = norm(max(abs(B(:,1:s)-Bstar(:,1:s)),[],1))^2;
        
        cur_kappa = LHS/RHS;
        kappa = min(kappa,cur_kappa);
    end
    kappa
    kappa_ls(ii) = kappa;
    
end