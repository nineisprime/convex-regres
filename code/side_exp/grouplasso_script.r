
# OUT:  
#   "gold" is smallest magnitude of correct
#   "lead" is the largest magnitude of incorr

runExp <- function(p, s, k, sigma, lambda){

	mu = rbind(matrix(abs(rnorm(s*k)),s,k), matrix(0,p-s,k));
	# mu = [non-zero; zero]
	
	epsilon = matrix(rnorm(p*k)*sigma, p,k);
	
	y = mu+epsilon;
	
	mu_hat = matrix(0,p,k);
	for (i in 1:p){
		mu_hat[i,] = y[i,] - l1Project(y[i,], lambda);
	}
	
	weights = apply(abs(mu_hat), 1, max);
	
	gold = min(weights[1:s]);
	lead = max(weights[(s+1):p]);
	
	return(list(gold=gold, lead=lead));
}

l1Project <- function(vec, lambda) {
	abs_ix = order(abs(vec),decreasing=TRUE);
	n = length(abs_ix);
	
	quantities = cumsum(abs( vec[abs_ix]) ) - 
				c( abs( vec[abs_ix[2:n]] ),0) *  1:n;
	
	if (quantities[n] < lambda) {
		return(vec);
	} else {	
		j_star = which (quantities >= lambda)[1];
	}
	
	
	if (j_star < n){
		vec[abs_ix[(j_star+1):n]] = 0;
	}
	current_sum = sum(abs(vec[abs_ix[1:j_star]]));
	
	vec[abs_ix[1:j_star]] = sign( vec[abs_ix[1:j_star]] ) * 
		(	abs(vec[abs_ix[1:j_star]]) -  (current_sum - lambda)/j_star );
		
	return(vec);
}	

k = 1000;
res = runExp(1000,3,k,0.5,k*0.05*log(k));
print(res)