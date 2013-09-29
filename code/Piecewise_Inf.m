function z = Piecewise_Inf(x,v)
n = length(x); [xa,ndx] = sort(abs(x),'descend');
b = [xa(2:n) 0] - (cumsum(xa)-v)./(1:n); z = zeros(1,n); 
if b(n) < 0
    pos = find(b<0,1); zinf = (sum(abs(x(ndx(1:pos))))-v)/pos;
    z(ndx(1:pos)) = sign(x(ndx(1:pos))).*zinf;
    z(ndx(pos+1:n)) = x(ndx(pos+1:n));
end
return