function R = calcrotmat(x)

a = angle(x(3) + 1i*x(4));

R = [cos(a), -sin(a); sin(a), cos(a)];

end