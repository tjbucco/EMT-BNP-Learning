function e = calcelipse(x, Sigma)

t = linspace(0,2*pi,100);
e0 = [cos(t) ; sin(t)];

R = calcrotmat(x);
[eigvec, eigval] = eig(round(R * Sigma * R.'));
eigvecsc = eigvec*sqrt(eigval);
e = kron(ones(1, 100), x(1:2)) + eigvecsc*e0;

end