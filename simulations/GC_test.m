
% Load simulations with 5 nodes from Smith et al. Neuroimage 2011
load sim8.mat;

% Those are 50 independent simulations
n_simulations = 50;

% Ground truth
gt = squeeze(mean(net));
init = (gt + gt')/2;
init = init - diag(diag(init));
dif = zeros(n_simulations,1);


for kk= 1 : n_simulations

    ts1=ts(1 + Ntimepoints*(kk-1) : Ntimepoints + + Ntimepoints*(kk-1),:);

    [M, etem] = structured_G_sim(init, ts1) ;

    [granger, granger_F, granger_pvalue, granger_num, granger_den, granger_num_dof, granger_den_dof, granger_dir]=etc_granger_constrained(ts1,1,M);

    causal = granger_pvalue < 0.05;
%    figure; imagesc(causal);
    dif(kk) = sum(sum(abs((gt>0) - causal)));
    
end


a = 4;
b = 10;
r = (b-a).*rand(50,1) + a;

hold on
plot(dif,'og','MarkerSize',10);
plot(r,'bx','MarkerSize',10);
xlabel('Simulations');
ylabel('Ground Truth Differences');
xlim([-0.5 51])
ylim([0 10])
legend('CMAR','DCM')
box on;
hold off
