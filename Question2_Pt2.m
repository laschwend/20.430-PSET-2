clear;
N = 2000; % simulate iterations
n = 6;
k = 1.5;
d = 0.1;
K = 2;

x1 = zeros(N,1); % assume the promoter is initially off arbitrarily
x2 = zeros(N,1); % column-vector for storing # mRNAs
x3 = zeros(N,1); % column-vector for storing # proteins
x3(1) = 10;
t = zeros(N,1); % column-vector for storing time
for iter = 2:N
    p_make_x1 = k / (1 + (x3(iter-1)/K)^n);
    p_degrade_x1 = d * x1(iter-1);
    p_make_x2 = k / (1 + (x1(iter-1)/K)^n);
    p_degrade_x2 = d * x2(iter-1);
    p_make_x3 = k / (1 + (x2(iter-1)/K)^n);
    p_degrade_x3 = d * x3(iter-1);
    reactions_total = p_make_x1 + p_degrade_x1 + p_make_x2 + p_degrade_x2 + p_make_x3 + p_degrade_x3;
    random_number = rand();

    x1(iter) = x1(iter-1);
    x2(iter) = x2(iter-1);
    x3(iter) = x3(iter-1);

    if random_number < p_make_x1/(reactions_total)
        x1(iter) = x1(iter-1) + 1;
    elseif random_number < (p_make_x1 + p_degrade_x1)/(reactions_total)
        x1(iter) = x1(iter-1) - 1;
    elseif random_number < (p_make_x1 + p_degrade_x1 + p_make_x2)/(reactions_total)
        x2(iter) = x2(iter-1) + 1;
    elseif random_number < (p_make_x1 + p_degrade_x1 + p_make_x2 + p_degrade_x2)/(reactions_total)
        x2(iter) = x2(iter-1) - 1;
    elseif random_number < (p_make_x1 + p_degrade_x1 + p_make_x2 + p_degrade_x2 + p_make_x3)/(reactions_total)
        x3(iter) = x3(iter-1) + 1;
    else
        x3(iter) = x3(iter-1) - 1;
    end
    t(iter) = t(iter-1) + exprnd(1/reactions_total);
end

figure; 
stairs(t,x1, 'k-', 'LineWidth', 1); xlim([t(1) t(end)]);
hold on;
stairs(t,x2, 'b-', 'LineWidth', 1); xlim([t(1) t(end)]);
hold on;
stairs(t,x3, 'r-', 'LineWidth', 1); xlim([t(1) t(end)]);
legend('$x_1$', '$x_2$', '$x_3$', 'Interpreter','latex')
xlabel('Time (seconds)');
ylabel('Number of Each Species')
hold off;
