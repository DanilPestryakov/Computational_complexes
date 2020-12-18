Aconst=[13 15;16 16;23 14];
A=[infsup(11,15) infsup(13,17);infsup(14,18) infsup(14,18);infsup(21,25) infsup(12, 16)];
%b = 9.7 11.2 13.4
x = [0.4; 0.3];
b=[infsup(7.2,12.2);infsup(9.2,13.2);infsup(10.9, 15.9)];
% solution with tolsolvty
% 3 x 2 matrix
inf_A = [11 13; 14 14; 21 12];
sup_A = [15 17; 18 18; 25 16];
inf_b = [7.2; 9.2; 10.9];
sup_b = [12.2; 13.2; 15.9];
[maxTol,argmaxTol,envs,ccode]=tolsolvty(inf_A,sup_A,inf_b,sup_b);
%Cminim = cond(inf_A);
b1=inf(A);
Cminim = condd(A, inf_A, sup_A, b1, 1);

b = 0.5 * (abs(sup_b)-abs(inf_b));
ive=sqrt(2)* Cminim * maxTol * norm(argmaxTol)/ norm(b);
disp('For the first task:');
disp('Ive_2D = ' + ive);
disp('Max tol = ' + maxTol);
disp('Cond = ' + Cminim);
cd('C:\Users\Данил\Documents\Intlab_V11\Intlab_V11\IntLinInc2D_1')
[V,P1,P2,P3,P4]=EqnTol2D(inf_A,sup_A,inf_b,sup_b);
cd('C:\Users\Данил\Documents\Intlab_V11\Intlab_V11')
rectangle('Position',[argmaxTol(1) argmaxTol(2) 0.001 0.001 ],'EdgeColor','r');
text(argmaxTol(1)+0.005,argmaxTol(2),'argmaxTol','FontSize',8);
title_str='3 x 2 task'; 
title(title_str);
xlabel('x_1');
ylabel('x_2');
title_str_name=strcat(title_str);
figure_name_out=strcat(title_str_name,'.png');
print('-dpng', '-r300', figure_name_out), pwd




% 2 x 3 matrix
%x = 0.4 0.3 0.2
%c = 14.6 13.6
c = [infsup(12.6, 16.6); infsup(11.1, 16.1)];
inf_c = [12.6; 11.1];
sup_c = [16.6; 16.1];
[maxTol2,argmaxTol2,envs,ccode] = tolsolvty(inf_A', sup_A', inf_c, sup_c);
b1 = inf(A');
Cminim2 = condd(A',inf_A', sup_A', b1, -1);
c = 0.5 * (abs(sup_c)-abs(inf_c));
ive2 = sqrt(3) * maxTol2 * norm(argmaxTol2) / norm(c) * Cminim2;
disp('For the second task:');
disp('Ive_2D = ' + ive2);
disp('Max tol = ' + maxTol2);
disp('Cond = ' + Cminim2);
cd('C:\Users\Данил\Documents\Intlab_V11\Intlab_V11\IntLinInc3D_1')
[V] = EqnTol3D(inf_A', sup_A', inf_c, sup_c, 1, 1);
cd('C:\Users\Данил\Documents\Intlab_V11\Intlab_V11')
title_str='2 x 3 task';
title(title_str);
xlabel('x_1');
ylabel('x_2');
zlabel('x_3');
title_str_name=strcat(title_str);
figure_name_out=strcat(title_str_name,'.png');
print('-dpng', '-r300', figure_name_out), pwd


Xsolv=[argmaxTol2-ive2/2,argmaxTol2+ive2/2];
A_c = [infsup(14,18) infsup(21,25); infsup(14,18) infsup(12, 16)];
cd('C:\Users\Данил\Documents\Intlab_V11\Intlab_V11\IntLinInc2D_1')
[V,P1,P2,P3,P4]=EqnTol2D(inf(A_c),sup(A_c),inf_c,sup_c);
cd('C:\Users\Данил\Documents\Intlab_V11\Intlab_V11')
rectangle('Position', [
Xsolv(2,1) Xsolv(3,1) Xsolv(2,2)-Xsolv(2,1) Xsolv(3,2)-Xsolv(3,1)]);
rectangle('Position',[argmaxTol2(2) argmaxTol2(3) 0.001 0.001 ],'EdgeColor','r');
title_str='Bar and tolerant set';
title(title_str);
xlabel('x_2');
ylabel('x_3');
xlim([0 1.2])
ylim([-0.33 0.67])

title_str_name=strcat(title_str);
figure_name_out=strcat(title_str_name,'.png');
print('-dpng', '-r300', figure_name_out), pwd
x = [0.4; 0.3; 0.2];
norm(x - argmaxTol2)