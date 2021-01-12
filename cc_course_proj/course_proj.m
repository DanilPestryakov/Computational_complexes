A = [3 5; -1 -3];
A(:, :, 2) = [4 6; 1 1];
b = [-3; -1];
b(:, :, 2) = [4; 2];
midb = [0.5; 0.5];
stop = 0;
for i=2.51:-0.001:2.49
  for j = 1.5:-0.01:0
    b = [midb(1, 1) - i; midb(2, 1) - j];
    b(:, :, 2) = [midb(1, 1) + i; midb(2, 1) + j];
    [result, opt, iter, epsilon] = subdiff(A, b);
    if iter < 20 || epsilon < 0.000001
      stop = 1;
      disp(result);
      disp(iter);
      disp(epsilon);
      disp(i);
      disp(j);
      break
    endif
  endfor
  if stop == 1
    break
  endif
endfor

for i=3.5:-0.1:3.4
  for j = 1.5:-0.1:0
    b = [midb(1, 1) - i; midb(2, 1) - j];
    b(:, :, 2) = [midb(1, 1) + i; midb(2, 1) + j];
    [result, opt, iter, epsilon] = subdiff(A, b);
    
    if (iter < 20 || epsilon < 0.000001)
      stop = 1;
      disp(result);
      disp(iter);
      disp(epsilon);
      disp(i);
      disp(j);
      
      break
    endif
  endfor
  if stop == 1
    break
  endif
endfor
b_new = b;
vec_dist(b, b_new);

new_rad_m = []
for i=0:1:10
  for j=0:1:10
    b(1, :, 1) = b_new(1, 1, 1) - i / 10;
    b(1, :, 2) = b_new(1, 1, 2) + j / 10;
    [result, opt, iter, epsilon] = subdiff(A, b);
    if (iter < 20 || epsilon < 0.000001)
      new_rad_m(i + 1, j + 1) = 1;
    else
      new_rad_m(i + 1, j + 1) = 0;
    endif
  endfor
endfor

b_new_1 = b_new;
b_new_1 = [-3; -1];
b_new_1(:, :, 2) = [3; 2];
vec_dist(b, b_new_1)
count = 0;
while 1
  b_new_1(1, 1, 2) = b_new_1(1, 1, 2) + 0.0001;
  [result, opt, iter, epsilon] = subdiff(A, b_new_1);
  if ~(iter < 20 || epsilon < 0.000001)
    b_new_1(1, 1, 2) = b_new_1(1, 1, 2) - 0.001;
    break
  endif 
  count = count + 1;
endwhile
b_new_1(1, 1, 2) = 3 + 0.0001;
[result, opt, iter, epsilon] = subdiff(A, b_new_1)


A = [3 5; -1 -3];
A(:, :, 2) = [4 6; 1 1];
b = [-3; -1];
b(:, :, 2) = [4; 2];
midA_2 = [0 -1];
[result, opt, iter, epsilon] = subdiff(A, b)
new_rad_m = []
for i=0:1:10
  for j=0:1:20
    A(2, 1, 1) = midA_2(1) - i / 10;
    A(2, 1, 2) = midA_2(1) + i / 10;
    A(2, 2, 1) = midA_2(2) - j / 10;
    A(2, 2, 2) = midA_2(2) + j / 10;   
    [result, opt, iter, epsilon] = subdiff(A, b);
    if (iter < 20 || epsilon < 0.000001)
      new_rad_m(i + 1, j + 1) = 1;
    else
      new_rad_m(i + 1, j + 1) = 0;
    endif
  endfor
endfor
imagesc(new_rad_m)
colorbar
xlabel('i')
ylabel('j')
title('Is have solution with rad((i - 1) / 10, (j - 1) / 10)')

dist_matr = [];
for i=0:1:10
  for j=0:1:20
    x_1_1 = midA_2(1) - i / 10;
    x_1_2 = midA_2(2) - j / 10;
    x_2_1 = -1;
    x_2_2 = -3;
    y_1_1 = midA_2(1) + i / 10;
    y_1_2 = midA_2(2) + j / 10;
    y_2_1 = 1;
    y_2_2 = 1;
    if new_rad_m(i + 1, j + 1) == 1
      dist_matr(i + 1, j + 1) = (sqrt((x_1_1 - x_2_1)^2 + (y_1_1 - y_2_1)^2) + sqrt((x_1_2 - x_2_2)^2 + (y_1_2 - y_2_2)^2)) / 2;
    else
      dist_matr(i + 1, j + 1) = 1000000;
    endif
  endfor
endfor

min_ij = 1000000;
i_m = 0;
i_j = 0;
for i=0:1:10
  for j=0:1:20
    if dist_matr(i + 1, j + 1) < min_ij
      min_ij = dist_matr(i + 1, j + 1);
      i_m = i;
      j_m = j;
    endif
  endfor
endfor

A_new = [3 5; -0.6 -3];
A_new(:, :, 2) = [4 6; 0.6 1];
b = [-3; -1];
b(:, :, 2) = [4; 2];
midA_2 = [0 -1];
[result, opt, iter, epsilon] = subdiff(A_new, b);

new_rad_A_m = []
for i=0:1:4
  for j=0:1:4
    A_new(2, 1, 1) = -0.6 - i / 10;
    A_new(2, 1, 2) = 0.6 + j / 10;
    [result, opt, iter, epsilon] = subdiff(A_new, b);
    if (iter < 20 || epsilon < 0.000001)
      new_rad_A_m(i + 1, j + 1) = 1;
    else
      new_rad_A_m(i + 1, j + 1) = 0;
    endif
  endfor
endfor

A_new = [3 5; -0.75 -3];
A_new(:, :, 2) = [4 6; 1 1];
b = [-3; -1];
b(:, :, 2) = [4; 2];
midA_2 = [0 -1];
[result, opt, iter, epsilon] = subdiff(A_new, b)
imagesc(new_rad_A_m)
colorbar
xlabel('i')
ylabel('j')
title('Deviation from [-0.6 0.6] with steps((i - 1) / 10, (j - 1) / 10)')