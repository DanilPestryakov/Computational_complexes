A = [2 -2; -2 2]
A(:, :, 2) = [4 1; 1 4]
b = [-2; -2]
b(:, :, 2) = [2; 2]
[result, opt, iter, epsilon] = subdiff(A, b)

A = [0.9 0.9 0.9; 0.9 0.9 0; 0.9 0 0]
A(:, :, 2) = [1.1 1.1 1.1; 1.1 1.1 0; 1.1 0 0]
b = [2.8; 2.2; 0.9]
b(:, :, 2) = [3.2; 1.8; 1.1]
[result, opt, iter, epsilon] = subdiff(A, b)

fid = fopen('matrix_n_phi_1.txt');
i = 1;
while ~(feof(fid))
  m = fgetl(fid);
  if ~isempty(m)
    A_128(i, :) = str2num(m);
    i = i + 1;
  endif
endwhile
fclose(fid);
s36 = [0.0 0.25 0.5 0.5 0.25 0.0 0.25 0.5 0.75 0.75 0.5 0.25 0.5 0.75 1.0 1.0 0.75 0.5 0.5 0.75 1.0 1.0 0.75 0.5 0.25 0.5 0.75 0.75 0.5 0.25 0.0 0.25 0.5 0.5 0.25 0.0];

A_128 = A_128(1:128, :);
deleting = [];
end_del = 1;
for i = 1:columns(A_128)
  count = 0;
  for j = 1:rows(A_128)
    if A_128(j, i) == 0
      count = count + 1;
    endif
    if count == rows(A_128)
      deleting(end_del) = i;
      end_del = end_del + 1;
    endif
  endfor
endfor
deleting
A_128(:, deleting) = [];
s36(:, deleting) = [];

norm_vec = []
for i = 1:rows(A_128)
  norm_vec(i) = norm(A_128(i, 1:18));
endfor
n = sort(norm_vec);

A_128_copy = A_128;
s36 = s36';
s36_copy = s36;
b_128 = A_128 * s36_copy;
b_128_copy = b_128;
dlmwrite('snew.txt', s36_copy)
dlmwrite('bnew.txt', b_128)

A_128 = A_128_copy;
b_128 = b_128_copy;
deleting = []
end_del = 1
for i = 1:rows(A_128)
  if norm_vec(i) < n(111)
    deleting(end_del) = i;
    end_del = end_del + 1;
  endif
endfor

A_128(deleting, :) = [];
b_128(deleting, :) = [];

A_18_1 = A_128;
A_18_2 = A_128;
rad_A = randn(18, 18);
b_18_1 = b_128;
b_18_2 = b_18_1;
rad_b = rand(18, 1);
A_18_3 = A_18_1 - rad_A;
A_18_3(:, :, 2) = A_18_2 + rad_A;
b_18_3 = b_18_1 - rad_b;
b_18_3(:, :, 2) = b_18_2 + rad_b;
[result, opt, iter, norm_diff] = subdiff(A_18_3, b_18_3)
x_opt_1 = result;
x_opt_2 = result;
x_opt_3 = result;
x_opt_4 = result;

A_128 = A_128_copy;
b_128 = b_128_copy;
deleting = []
end_del = 1
for i = 1:rows(A_128)
  if ~(norm_vec(i) < n(111) && norm_vec(i) > n(92))
    deleting(end_del) = i;
    end_del = end_del + 1;
  endif
endfor

A_128(deleting, :) = [];
b_128(deleting, :) = [];

A_18_1 = A_128;
A_18_2 = A_128;
rad_A = randn(18, 18);
b_18_1 = b_128;
b_18_2 = b_18_1;
rad_b = rand(18, 1);
A_18_3 = A_18_1 - rad_A;
A_18_3(:, :, 2) = A_18_2 + rad_A;
b_18_3 = b_18_1 - rad_b;
b_18_3(:, :, 2) = b_18_2 + rad_b;
[result, opt, iter, norm_diff] = subdiff(A_18_3, b_18_3)
x_opt_2 = result


A_128 = A_128_copy;
b_128 = b_128_copy;
deleting = []
end_del = 1
for i = 1:rows(A_128)
  if ~(norm_vec(i) < n(93) && norm_vec(i) > n(74))
    deleting(end_del) = i;
    end_del = end_del + 1;
  endif
endfor

A_128(deleting, :) = [];
b_128(deleting, :) = [];

A_18_1 = A_128;
A_18_2 = A_128;
rad_A = randn(18, 18);
b_18_1 = b_128;
b_18_2 = b_18_1;
rad_b = rand(18, 1);
A_18_3 = A_18_1 - rad_A;
A_18_3(:, :, 2) = A_18_2 + rad_A;
b_18_3 = b_18_1 - rad_b;
b_18_3(:, :, 2) = b_18_2 + rad_b;
[result, opt, iter, norm_diff] = subdiff(A_18_3, b_18_3)
x_opt_3 = result


A_128 = A_128_copy;
b_128 = b_128_copy;
deleting = []
end_del = 1
for i = 1:rows(A_128)
  if (norm_vec(i) > n(18))
    deleting(end_del) = i;
    end_del = end_del + 1;
  endif
endfor

A_128(deleting, :) = [];
b_128(deleting, :) = [];

A_18_1 = A_128;
A_18_2 = A_128;
rad_A = randn(18, 18);
b_18_1 = b_128;
b_18_2 = b_18_1;
rad_b = rand(18, 1);
A_18_3 = A_18_1 - rad_A;
A_18_3(:, :, 2) = A_18_2 + rad_A;
b_18_3 = b_18_1 - rad_b;
b_18_3(:, :, 2) = b_18_2 + rad_b;
[result, opt, iter, norm_diff] = subdiff(A_18_3, b_18_3)
x_opt_4 = result

x_opt_int_inf = [];
x_opt_int_sup = [];
x_opt_int_mid = [];
for i = 1:18
  x_opt_int_inf(i) = max([x_opt_1(i) x_opt_2(i) x_opt_3(i) x_opt_4(i)]);
  x_opt_int_sup(i) = min([x_opt_1(i + 18) x_opt_2(i + 18) x_opt_3(i + 18) x_opt_4(i + 18)]);
  x_opt_int_mid(i) = (x_opt_int_inf(i) + x_opt_int_sup(i)) / 2;
endfor
plot(x_opt_int_inf, 'r')
hold on
plot(x_opt_int_sup, 'b')
title('18x18 solve')
xlabel('i')
ylabel('x[i]')
legend('inf', 'sup')



b_mid = []
for i = 1:18
  b_mid(i) = 0;
  for j = 1:18
    t = times([A_18_3(i, j, 1) A_18_3(i, j, 2)], [x_opt_int_inf(j) x_opt_int_sup(j)]);
    b_mid(i) = b_mid(i) + (t(1, 1) + t(1, 2)) / 2;
  endfor
endfor
plot(b_18_3(:, :, 1), 'r')
hold on
plot(b_18_3(:, :, 2), 'b')
hold on
plot(b_mid, 'g--')
xlabel('i')
ylabel('b[i]')
title('18x18 vector b')
legend('b_{inf}', 'b_{sup}', 'midAx')