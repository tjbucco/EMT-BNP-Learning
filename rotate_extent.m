function D = rotate_extent(x, D)

for i = 1:size(x, 2)
a = angle(x(3, i) + 1i*x(4, i));

R = [cos(a), -sin(a); sin(a), cos(a)];

D(:,:,i) = R * D(:,:,i) * R.';
end
end