clear;
N = 20000; % simulate iterations

%Define rates
k_1 = 0.05; 
d_1 = 0.11;
k_2 = 5;
d_2 = 1;
k_3 = 20;
d_3 = 0.1;

P = zeros(N,1); % assume the promoter is initially off arbitrarily
mRNA = zeros(N,1); % column-vector for storing # mRNAs
protein = zeros(N,1); % column-vector for storing # proteins
t = zeros(N,1); % column-vector for storing time
for iter = 2:N
    if P(iter-1) == 0
        p_promoter = k_1;
    else
        p_promoter = d_1;
    end
    p_make_mRNA = k_2 * P(iter-1);
    p_degrade_mRNA = d_2 * mRNA(iter-1);
    p_make_protein = k_3 * mRNA(iter-1);
    p_degrade_protein = d_3 * protein(iter-1);

    reactions_total = p_promoter + p_make_mRNA + p_degrade_mRNA + p_make_protein + p_degrade_protein;
    random_number = rand();
    if random_number < p_promoter/(reactions_total)
        if P(iter-1) == 0
            P(iter) = 1;
            mRNA(iter) = mRNA(iter-1);
            protein(iter) = protein(iter-1);
        else
            P(iter) = 0;
            mRNA(iter) = mRNA(iter-1);
            protein(iter) = protein(iter-1);
        end
    elseif random_number < (p_promoter + p_make_mRNA)/(reactions_total)
        mRNA(iter) = mRNA(iter-1) + 1;
        P(iter) = P(iter-1);
        protein(iter) = protein(iter-1);
    elseif random_number < (p_promoter + p_make_mRNA + p_degrade_mRNA)/(reactions_total)
        mRNA(iter) = mRNA(iter-1) - 1;
        P(iter) = P(iter-1);
        protein(iter) = protein(iter-1);
    elseif random_number < (p_promoter + p_make_mRNA + p_degrade_mRNA + p_make_protein)/(reactions_total)
        protein(iter) = protein(iter-1) + 1;
        mRNA(iter) = mRNA(iter-1);
        P(iter) = P(iter-1);
    else
        protein(iter) = protein(iter-1) - 1;
        mRNA(iter) = mRNA(iter-1);
        P(iter) = P(iter-1);
    end
    t(iter) = t(iter-1) + exprnd(1/reactions_total);
end

figure;
set(gcf, 'Units', 'inches', 'Position', [1, 1, 10, 4]);
subplot(1,3,1);
stairs(t,P, 'k-', 'LineWidth', 1); xlim([t(1) t(end)]);
ylabel('Promoter state'); xlabel('time (seconds)');

subplot(1,3,2);
stairs(t,mRNA, 'b-', 'LineWidth', 1); xlim([t(1) t(end)]);
ylabel('mRNA number'); xlabel('time (seconds)');

subplot(1,3,3);
stairs(t,protein, 'r-', 'LineWidth', 1); xlim([t(1) t(end)]);
ylabel('Protein number'); xlabel('time (seconds)');
