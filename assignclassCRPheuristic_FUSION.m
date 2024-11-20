function [D, Q, Cnew] = assignclassCRPheuristic_FUSION(newcust_total, D,Q,C,hyper)

newcust_count = 0;
prob_existing_table_ell = NaN(1, size(C, 2) + newcust_total);
Cnew = zeros(newcust_total, 1);
Dnew = NaN(2, 2, newcust_total);
Sigmavnew = NaN(4, 4, newcust_total);

while newcust_count < newcust_total

    total_customers = size(C,2) + newcust_count;
    for l = 1:(max([C.', Cnew.'].'))
        customer_count = sum([C.', Cnew.'].' == l);
        prob_existing_table_ell(l) = customer_count/(hyper.alphadp + total_customers);
    end
    if hyper.alphadp == inf
        prob_existing_table_ell(1:l) = 0; 
        prob_new_table = 1;
    else
        prob_new_table = hyper.alphadp/(hyper.alphadp + total_customers);
    end
    newcust_count = newcust_count + 1;

    Cnew(newcust_count) = find(mnrnd(1, [prob_existing_table_ell(1:l), prob_new_table]./(sum([prob_existing_table_ell(1:l), prob_new_table]))));

    % Generate extent matrices
    lsqinv = gamrnd(hyper.al0, 1./hyper.bl0, 2, 1);
    lsq = 1./lsqinv;
    Dnew(:,:, newcust_count) = diag(lsq);

    % Generate driving noise covariance matrices
    ssqinv = gamrnd(hyper.av0, 1./hyper.bv0, 2, 1);
    ssq = 1./ssqinv;
    Sigmavnew(:, :, newcust_count) = diag(ssq([1, 1, 2, 2]));

end

%C = [C, Cnew];
D = cat(3, D, Dnew);
Q = cat(3, Q, Sigmavnew);
end
