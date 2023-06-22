%% Exercíco programa 3 - Exercicio 1
% ----------------------------------------------------------------
% Data: Maio/2022
% Autor: Bruno Fernandes
% ----------------------------------------------------------------

%Carregar a série "pib.csv"

rng(123)

% Construindo o vetor de datas trimestrais de 1999 a 2021
count = 1;
for i = 1999:2021
    for j = 1:3:12
        data(count) = datetime(i, j, 01);
        count = count + 1;
    end
end

data.Format = 'mmm-yyyy';

% Log do PIB dessazonalizado de jan/1999-out/2021
PIB = PIB(13:104);

figure(7)
plot(data, PIB);
title('Log do PIB dessazonalizado de jan/1999-out/2021')

%Estabelecendo as priors dos parâmetros
phi_pr = [0.7; 0]; % média de phi
Vp_pr = diag([0.1 0.05]); % variância de phi

delta_pr = 0.005; % média de delta
Vd_pr = 0.01; %variância de delta

nuv_pr = 10; % nuv_pr/2 = alfa da distribuição inversa gamma da variância de V
nuw_pr = 10; % sv_pr/2 = beta da distribuição inversa gamma de variância de V

sv_pr = 0.1; % nuw_pr/2 = alfa da distribuição inversa gamma de variância de W
sw_pr = 0.1; % sw_pr/2 = beta da distribuição inversa gamma de variância de W

eta_pr = [0; 0]; % média de eta
p_pr = 0.01*eye(2); % variância de eta

%% a)

% Plot das priors
grid_norm = -2:0.01:2;
grid_IG = 0.001:0.001:0.1;

figure(1)
subplot(2,4,1)
    hold on;
    plot(grid_norm, normpdf(grid_norm, phi_pr(1), sqrt(Vp_pr(1, 1))),'LineWidth', 2)
    xlabel('Phi1')
subplot(2,4,2)
    hold on;
    plot(grid_norm, normpdf(grid_norm, phi_pr(2), sqrt(Vp_pr(2, 2))), 'LineWidth', 2)
    xlabel('Phi2')
subplot(2,4,3)
    hold on;
    plot(grid_norm, normpdf(grid_norm, delta_pr, sqrt(Vd_pr)), 'LineWidth', 2)
    xlabel('Delta')
subplot(2,4,4)
    hold on;
    plot(grid_IG, gampdf(1./grid_IG, nuv_pr/2, 2/sv_pr), 'LineWidth', 2)
    xlabel('Sigma^2_v')
subplot(2,4,5)
    hold on;
    plot(grid_IG, gampdf(1./grid_IG, nuw_pr/2, 2/sw_pr), 'LineWidth', 2)
    xlabel('Sigma^2_w')
subplot(2,4,6)
    hold on;
    plot(grid_norm, normpdf(grid_norm, eta_pr(1), sqrt(p_pr(1, 1))), 'LineWidth', 2)
    xlabel('Eta_0')
subplot(2,4,7)
    hold on;
    plot(grid_norm, normpdf(grid_norm, eta_pr(2), sqrt(p_pr(2, 2))), 'LineWidth', 2)
    xlabel('Eta_1')
    
%% d)

draws = 7000;
burnin = 2000;
T = size(PIB,1);

% Matrizes de armanezamento
phi = nan(draws, 2); 
sigw = nan(draws, 1); %variância de W

delta = nan(draws, 1);
sigv = nan(draws, 1); %variância de V

eta = nan(T-1, 2, draws);

% Valores iniciais (dados pela média de cada parâmetro)
phi(1,:) = phi_pr';
sigw(1) = (sw_pr/2)/((nuw_pr/2) - 1);
delta(1) = delta_pr;
sigv(1) = (sv_pr/2)/((nuv_pr/2) - 1);

% Matrizes iniciais do espaço estado
F = [phi(1, :); 1 0];
G = [1; 0];
Q = sigw(1);
H0 = delta(1);
H1 = [1 -1];
R = sigv(1);

% valores posteriors dos hiperparâmetros nuv_ps e nuw_ps
nuv_ps = nuv_pr + T;
nuw_ps = nuw_pr + T;

% Tirando a primeira diferença do log do Pib

y = PIB(2:end) - PIB(1:end-1);

% Algorítimo de sorteio da posterior com modelo de espaço estado linear
for i = 2:draws
    eta(:,:, i) = carter_kohn_sampler(F, G, Q, H0, H1, R, y, eta_pr, p_pr);
    
    %Atualização dos hiperparâmetros para sorteio de novos parâmetros da posterior
    Y = eta(2:end, 1, i);
    X = eta(1:end-1, :, i);
    phi_chap = inv(X'*X)*(X'*Y);
    
    %Atualização de sigw
    sw_ps = sw_pr + (Y - X*phi(i-1,:)')'*(Y - X*phi(i-1,:)');
    sigw(i) = 1./gamrnd(nuw_ps/2, 2/sw_ps);
    Q = sigw(i);
    
    %Atualização de phi
    Vp_ps = inv(inv(Vp_pr) + (1/sigw(i))*(X'*X));
    phi_ps = Vp_ps*(inv(Vp_pr)*phi_pr + (1/sigw(i))*(X'*X)*phi_chap); 
    phi(i,:) = mvnrnd(phi_ps, Vp_ps, 1);
    F = [phi(i, :); 1 0];

    mu = PIB(2:end) - eta(:, 1, i);
    Y = mu(2:end) - mu(1:end-1);
    X = ones(T-2, 1);
    delta_chap = inv(X'*X)*(X'*Y);
    
    %Atualização de sigv
    sv_ps = sv_pr + (Y - X*delta(i-1))'*(Y - X*delta(i-1));
    sigv(i) = 1./gamrnd(nuv_ps/2, 2/sv_ps);
    R = sigv(i);

    %Atualização de delta
    Vd_ps = inv(inv(Vd_pr) + (1/sigv(i))*(X'*X));
    delta_ps = Vd_ps*(inv(Vd_pr)*delta_pr + (1/sigv(i))*(X'*X)*delta_chap); 
    delta(i) = mvnrnd(delta_ps, Vd_ps, 1);
    H0 = delta(i);

end

% Burn-ins
phi = phi(burnin+1:end,:); 
sigw = sigw(burnin+1:end);
delta = delta(burnin+1:end);
sigv = sigv(burnin+1:end);
eta = eta(:,:,burnin+1:end);

theta = [phi sigw delta sigv];
mean_ps = mean(theta)
std_ps = std(theta)

%% e)

% Gráfico Histograma
var_name = {'phi_1' 'phi_2' 'sigma^2_w' 'delta' 'sigma^2_v'};
k = 5;

figure(2)
set(gca, 'Color', 'w');
for i = 1:k
    subplot(2, 3, i)
    histogram(theta(:, i), 100, 'FaceColor', '#0072BD')
    hold on
    ylimits = ylim;
    plot([mean_ps(i) mean_ps(i)], [ylimits(1) ylimits(2)], 'LineWidth', 3)
    title(var_name{i})
end


% Gráfico autocorrelações
% Função autocorr proveniente do pacote Econometrics dentro da licença da GV 
n_acf = 50;
window = 50;
n_draws = draws - burnin;

figure(3)
set(gca, 'Color', 'none');
for i = 1:k % for each variable
    [acf,lag,bounds] = autocorr(theta(:,i),n_acf); 
    movavg = movmean(theta(:,i),window); 
    s(i) = subplot(2,k,i);
    plot(theta(:,i))
    hold on;
    subplot(2,k,i), plot(1:n_draws,movavg,'Color','r','LineWidth', 2);
    xlim([0, n_draws + 1])
    title(var_name{i})

    subplot(2,k,k+i),stem(acf,'Color','r');
    subplot(2,k,k+i),line([min(lag) max(lag)], [min(bounds) min(bounds)],'Color','b','LineWidth',1);
    subplot(2,k,k+i),line([min(lag) max(lag)], [max(bounds) max(bounds)],'Color','b','LineWidth',1);
    xlim([0, n_acf + 1])
    title('ACF')
end
legend(s(3),'Sorteios','Média Móvel (50 sorteios)')

%% e)

% Squeeze de mu e eta
mu = squeeze(PIB(2:end) - eta(:, 1,:));
eta = squeeze(eta(:, 1,:));

% Gráfico produto vs. produto Potencial
figure(4)
set(gca, 'Color', 'w');
hold on
plot(data(2:end), PIB(2:end), 'k', 'LineWidth', 2)
plot(data(2:end), mean(mu, 2), 'b')
shade_interval(data(2:end), quantile(mu, 0.16, 2), 'b--', data(2:end), quantile(mu, 0.84, 2), '--b', 'FillType', [1 2; 2 1], 'FillAlpha', 0.1)
legend('Produto','Produto Potencial')
grid on
grid minor

%% f)

% Gráfico hiato do produto
figure(5)
set(gca, 'Color', 'w');
hold on
plot(data(2:end), zeros(T-1), 'k-.', 'LineWidth', 1.5)
plot(data(2:end), mean(eta, 2), 'b-', 'LineWidth', 1.5)
shade_interval(data(2:end), quantile(eta, 0.16, 2), 'r:', data(2:end), quantile(eta, 0.84, 2), 'r:', 'FillType', [1 2; 2 1], 'FillAlpha', 0.1)
grid on
grid minor


