x1s = (0:1e2)/1e2;
x2s = (0:1e2)/1e2;

f = zeros(length(x1s),length(x2s));

for ii = 1:length(x1s)
    for jj=1:length(x2s)
        x1 = x1s(ii);
        x2 = x2s(jj);
        
        if (x1 < x2)
            f(ii,jj) = -2*(x1/x2) + 1;
        else 
            f(ii,jj) = x1*(2/(1-x2)) - 2/(1-x2) + 1;
        end
    end
end

