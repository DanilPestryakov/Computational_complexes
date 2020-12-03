% bar for solution
X = intval([infsup(-5, 5), infsup(-5, 5)]);
% rasstrigin function
[Z, WorkList, diams] = globopt0(X);
iter = 1:1:length(diams);
plot(iter, diams);
hold on;
xlim([0, length(diams)]);
xlabel('Iterations');
ylabel('Diameter of area');
title('Constriction of area');
path = 'C:\Octave\Intlab_V11\Intlab_V11';
full_title = 'Rasstrigin';
saveas(gcf, fullfile(path, char(full_title)), 'png'); 

for i = 1:30
    disp(WorkList(i).Box);
    s = ['f(y) = ', num2str(WorkList(i).Estim)];
    disp(s);
end

%McCormick function
X = intval([infsup(-1.5, 4), infsup(-3, 4)]);
[Z, WorkList, diams] = globopt0(X);
solution = -1.9133;

answer = [];
for i = 1 : length(WorkList)
    answer(i) = WorkList(i).Estim;
end

for i = 1:length(answer)
    diff(i) = abs(answer(i) - solution);
end

iter = 1:1:length(answer);
plot(iter, diff);

semilogx(iter, diff);
hold on;
xlim([0, length(answer)]);
xlabel('Iterations');
ylabel('Absolute difference');
title('Convergence of method');
path = 'C:\Octave\Intlab_V11\Intlab_V11';
full_title = 'McCormick function logx';
saveas(gcf, fullfile(path, char(full_title)), 'png'); 
WorkList(length(WorkList)).Box

% Bar's center
centers_x= [];
centers_y = [];
for i = 1 : length(WorkList)
    centers_x(i) = WorkList(i).Box(1).mid;
    centers_y(i) = WorkList(i).Box(2).mid;
end
plot(centers_x(900:1001), centers_y(900:1001)); 
xlabel('Center X');
ylabel('Center Y');
title('Trajectory of bar center');
full_title = 'Trajectory_center';
saveas(gcf, fullfile(path, char(full_title)), 'png'); 