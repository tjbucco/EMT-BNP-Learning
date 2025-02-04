function Xsf = csmcsampleauxpfas(Zf, Xcf, Ts, Te, D, F, Sigmav, Sigmax1, R)

Z = Zf(:, :, Ts:Te);
Xc = Xcf(:, Ts:Te);

maxno = size(Z, 3);
no = sum(~isnan(Xc(1, :)), 2);
maxnm = size(Z, 2);

M = NaN(1, maxno);
for t = 1:no
    if all(isnan(Z(1, :, t)))
        M(1, t) = 0;
    else
        M(1, t) = find(~isnan(Z(1, :, t)), 1, 'last');
    end
end

maxlsq = max(diag(D));
sv1sq = Sigmav(1, 1);
sv2sq = Sigmav(3, 3);

X = NaN(4, R, no);
W = NaN(1, R, no);
A = NaN(1, R, no - 1);

X(:, R, :) = Xc(:, 1:no);


Znm = NaN(2, maxno);
Znm(:, :) = nanmean(Z(:, :, :), 2);

% use weight exponent as unormalized weights can be small
wexpo = NaN(1, R);
wasexpo = NaN(1, R);

% first timestep
if M(1, 1)
    Sigmax1prime = blkdiag(eye(2)*max(diag(D))/M(1, 1), Sigmax1(3:4, 3:4));
    
    %X(:, 1:(R - 1), 1) = mvnrnd([Znm(:, 1).', 0, 0], Sigmax1prime, R - 1).';
    X(:, 1:(R - 1), 1) = mvnrnd_fast(repmat([Znm(:, 1).', 0, 0],R-1,1), diag(Sigmax1prime).').';

    for r = 1:R
        Rkt = calcrotmat(X(:, r, 1));
        Ynksqsum = nansum((Rkt.'*(Z(:, :, 1) - X(1:2, r*ones(1, maxnm), 1))).^2, 2);
        
        wexpo(r) = 1/2*(- 1./diag(Sigmax1).'*X(:, r, 1).^2 - 1./diag(D).'*Ynksqsum + 1./diag(Sigmax1prime).'*(X(:, r, 1) - [Znm(:, 1); 0; 0]).^2);
    end
else
    X(:, 1:(R - 1), 1) = mvnrnd([0, 0, 0, 0], Sigmax1, R - 1).';
    
    wexpo = ones(1, R);
end

wexpo = wexpo - max(wexpo);
w = exp(wexpo);
W(1, :, 1) = w/sum(w);

for t = 2:no
    [~, Atemp] = histc(rand(1, R - 1), [0, cumsum(W(1, :, t - 1))]);
    A(1, 1:(R - 1), t - 1) = Atemp;
    
    for r = 1:(R - 1)
        if ~M(1, t)
            X(:, r, t) = F*X(:, A(1, r, t - 1), t - 1) + mvnrnd_fast([0, 0, 0, 0], diag(Sigmav).').';
        else
            rho = sv1sq/(sv1sq + maxlsq/M(1, t));
            Xkp1 = F*X(:, A(1, r, t - 1), t - 1);
            
            uq = Xkp1 + rho*[sum(Z(:, 1:M(1, t), t), 2)/M(1, t) - Xkp1(1:2); 0; 0];
            Sigmaq = Sigmav - sv1sq*rho*diag([1, 1, 0, 0]);
            
            X(:, r, t) = uq + mvnrnd_fast([0, 0, 0, 0], diag(Sigmaq).').';

        end
    end
    
    for r = 1:R
        v = X(:, R, t) - F*X(:, r, t - 1);
        
        wasexpo(r) = wexpo(r) - 1/2*(1./diag(Sigmav).'*v.^2);
    end
    
    wasexpo = wasexpo - max(wasexpo);
    was = exp(wasexpo);
    Was = was/sum(was);
    
    [~, Atemp] = histc(rand, [0, cumsum(Was)]);
    A(1, R, t - 1) = Atemp;
    
    wexpo = NaN(1, R);
    
    for r = 1:R
        rho = sv1sq/(sv1sq + maxlsq/M(1, t));
        Xkp1 = F*X(:, A(1, r, t - 1), t - 1);
        
        uq = Xkp1 + rho*[sum(Z(:, 1:M(1, t), t), 2)/M(1, t) - Xkp1(1:2); 0; 0];

        Rkt = calcrotmat(X(:, r, t));
        Yktnsqsum = nansum((Rkt.'*(Z(:, :, t) - X(1:2, r*ones(1, maxnm), t))).^2, 2);
        
        if ~M(1, t)
            wexpo(r) = -1/2*1./diag(D).'*Yktnsqsum;
        else
            wexpo(r) = 1/2*(2*log(1 - rho) - 1./[sv1sq*[1, 1], sv2sq*[1, 1]]*(X(:, r, t) - Xkp1).^2 + 1./[sv1sq*(1 - rho)*[1, 1], sv2sq*[1, 1]]*(X(:, r, t) - uq).^2 - 1./diag(D).'*Yktnsqsum);
        end
    end
    
    wexpo = wexpo - max(wexpo);
    w = exp(wexpo);
    W(1, :, t) = w/sum(w);
end

Bs = NaN(1, no);
Xs = NaN(4, maxno);

[~, Bstemp] = histc(rand, [0, cumsum(W(1, :, no))]);
Bs(no) = Bstemp;
Xs(:, no) = X(:, Bs(no), no);

for t = fliplr(1:(no - 1))
    Bs(t) = A(1, Bs(t + 1), t);
    Xs(:, t) = X(:, Bs(t), t);
end

Xsf = NaN(4, size(Xcf, 2));
Xsf(:, Ts:Te) = Xs;



end