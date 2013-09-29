#IN:  "X"  (k--by--p) matrix 
#		"Y"  (k) vector,		"lambda"
#		
#OUT: "B"  (p--by--k) matrix
#
# solves:
# 1/2 | Y - diag(X*B) |^2 + lambda | B |_{1,inf}

multiLasso <- function(X, Y, lambda){

	k = nrow(X);  p=ncol(X);  max_iter = 200;
	B = matrix(0, p, k);
	
	prev_obj = Inf;
	stepsize = 0.5;
	armijo_sigma = 0.2;
	armijo_beta = 0.9;
	
	for (it in 1:max_iter){
		
		# loss_gradient
		loss_gradient = - diag( as.vector( Y - diag(X %*% B )) )%*% X;
		loss_gradient = t(loss_gradient); # p--by--k after tranpose
		
		# penalty_gradient
		#max_ixs = apply(abs(B), 1, which.max);
		#max_ixs2 = cbind(1:p, max_ixs);
		#pen_gradient = matrix(0, p, k);
		#pen_gradient[max_ixs2] = lambda*sign(B[max_ixs2]) * matrix(1,p,1);
		
		new_B = matrix(0, p, k);
		new_obj = 0;
		while(TRUE){
			tmp = B - stepsize*loss_gradient;
			for (ii in 1:p){
				new_B[ii,] = tmp[ii,] - l1Project(tmp[ii,], lambda);
			}
			
			new_obj = 0.5 * sum( (Y - diag(X %*% new_B)) ^ 2) + lambda * sum(apply(abs(new_B), 1, max));
			
			if ( new_obj < prev_obj + 
			          sum(diag( t(loss_gradient) %*% (new_B - B) )) + 
					   1/(2*stepsize) * sum(apply((new_B - B)^2,1,sum)) ){
				break;
			} else {
				stepsize = stepsize*armijo_beta;
			}
			
			if (stepsize < 1e-5){
				print(c(prev_obj, new_obj));
				print("Line search does not terminate");
				return(0);
			}
		} #end-while
		
		B = new_B;		
		print(c("iter:", it, "obj:", new_obj, "stepsize:", stepsize)); 
		flush.console();
		if ( abs(prev_obj - new_obj) < 1e-5*abs(prev_obj) ){
			break;
		}
		prev_obj = new_obj;
		
	} # end-for
	
	return(B);
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

# test case
k = 1000; p=50; s=3; noise = 0;
X = matrix(rnorm(k*p), k, p);
B_star = rbind(matrix(rnorm(s*k),s,k), matrix(0,(p-s),k));
Y = diag(X %*% B_star) + matrix(rnorm(k),k,1)*noise;

lambda = 4;
B = multiLasso(X,Y,lambda);
weights = apply(abs(B),1,max);
min(weights[1:s])
max(weights[(s+1):p])