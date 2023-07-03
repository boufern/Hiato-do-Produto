%% Output Gap
% -------------------------------------------------------------------------------------
% 
% Purpose: Bayesian estimation of a State-Space model that is able to
% estimate the economy's actual business cycle and potential output.
%
% Target dynamic model with constant parameters:
% y(t) = mu(t) + eta(t)
% mu(t) = delta + mu(t-1) + w(t)                   , w(t)~N(0,Q) iid
% eta(t) = phi1*eta(t-1) + phi2*eta(t-2) + v(t)    , v(t)~N(0,R) iid
% 
% where y(t) -- log of output
%       mu(t) -- potential outup
%       eta(t) -- output gap    
%
% State-space model with constant parameters used for estimation
% y(t) = H0 + H1*eta(t) + w(t),  w(t)~N(0,Q) iid
% eta(t) = F*eta(t-1) + G*v(t) ,  v(t)~N(0,R) iid
% 
% Inputs:   F  -- 2 x 2 matrix
%           G  -- 2 x 1 matrix
%           Q  -- 1 x 1 covariance matrix
%           H0 -- 1 x 1 vector of constants
%           H1 -- 2 x 1 matrix
%           R  -- 1 x 1 covariance matrix
%           y  -- T x 1 data matrix
% Outputs:  mu -- (T-1) x (draws-burnin) matrix
%           eta -- (T-1) x (draws-burnin) matrix
%
% Priors and        phi1 ~ N(0.7,0.1)
% initial values:   phi2 ~ N(0,0.05)
%                   delta ~ N(0.005,0.01)
%                   sigma_w ~ IG(5,0.05)
%                   sigma_v ~ IG(5,0.05)
%                   eta(0) ~ N([0 0], 0.01*I(2))

% -------------------------------------------------------------------------------------
% Date: 22/06/2023
% Author: Bruno Fernandes
% If you find any error, please contact bruno.fernandes03@outlook.com
% -------------------------------------------------------------------------------------

PIB = readmatrix("~/Files/Github/Hiato do Produto/PIB.csv");
PIB = PIB(2:end,3);
[T,~] = size(PIB);

rng(123)

% Building the vector of dates 
count = 1;
for i = 1996:floor(1996+T/4)
    for j = 1:3:12
        date(count) = datetime(i, j, 01);
        count = count + 1;
    end
end

formatOut = 'mmm-yyyy';
datestr(date,formatOut);

date = date(1:T);

% Setting up the parameter priors
phi_pr = [0.7; 0]; % mean of phi
Vp_pr = diag([0.1 0.05]); % variance of phi

delta_pr = 0.005; % mean of delta
Vd_pr = 0.01; % variance of delta

nuv_pr = 10; % nuv_pr/2 =  Inverse Gamma distribution alpha of the variance of V
nuw_pr = 10; % sv_pr/2 = Inverse Gamma distribution beta of the variância of V

sv_pr = 0.01*nuv_pr; % nuw_pr/2 = Inverse Gamma distribution alpha of the variance of W
sw_pr = 0.01*nuw_pr; % sw_pr/2 = Inverse Gamma distribution beta of the variância of W

eta_pr = [0; 0]; % mean eta
p_pr = 0.01*eye(2); % mean of eta

%Setting up the paramerters of the estimation
draws = 7000;
burnin = 2000;

% Storage matrixes 
phi = nan(draws, 2); 
sigw = nan(draws, 1); % variance of W

delta = nan(draws, 1);
sigv = nan(draws, 1); % variance of V

eta = nan(T-1, 2, draws);

% Initial Values (given by the mean of each parameter)
phi(1,:) = phi_pr';
sigw(1) = (sw_pr/2)/((nuw_pr/2) - 1);
delta(1) = delta_pr;
sigv(1) = (sv_pr/2)/((nuv_pr/2) - 1);

% Initial matrixes of the state-space model
F = [phi(1, :); 1 0];
G = [1; 0];
Q = sigw(1);
H0 = delta(1);
H1 = [1 -1];
R = sigv(1);

% The hyperparameters posterior values nuv_ps and nuw_ps
nuv_ps = nuv_pr + T;
nuw_ps = nuw_pr + T;

% Log of Pib first differences

y = PIB(2:end) - PIB(1:end-1);

% Drawing from the state space posterior
for i = 2:draws
    eta(:,:, i) = carter_kohn_sampler(F, G, Q, H0, H1, R, y, eta_pr, p_pr);
    
    % Hyperparameters updating in order to draw new parameters from the posterior
    Y = eta(2:end, 1, i);
    X = eta(1:end-1, :, i);
    phi_chap = inv(X'*X)*(X'*Y);
    
    % Updating sigw
    sw_ps = sw_pr + (Y - X*phi(i-1,:)')'*(Y - X*phi(i-1,:)');
    sigw(i) = 1./gamrnd(nuw_ps/2, 2/sw_ps);
    Q = sigw(i);
    
    % Updating phi
    Vp_ps = inv(inv(Vp_pr) + (1/sigw(i))*(X'*X));
    phi_ps = Vp_ps*(inv(Vp_pr)*phi_pr + (1/sigw(i))*(X'*X)*phi_chap); 
    phi(i,:) = mvnrnd(phi_ps, Vp_ps, 1);
    F = [phi(i, :); 1 0];

    mu = PIB(2:end) - eta(:, 1, i);
    Y = mu(2:end) - mu(1:end-1);
    X = ones(T-2, 1);
    delta_chap = inv(X'*X)*(X'*Y);
    
    % Updating sigv
    sv_ps = sv_pr + (Y - X*delta(i-1))'*(Y - X*delta(i-1));
    sigv(i) = 1./gamrnd(nuv_ps/2, 2/sv_ps);
    R = sigv(i);

    % Updating delta
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


% Autocorrelation Graphic
n_acf = 50;
window = 50;
n_draws = draws - burnin;

figure(1)
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

% Squeeze mu and eta
mu = squeeze(PIB(2:end) - eta(:, 1,:));
eta = squeeze(eta(:, 1,:));

% Actual output vs. the potencial one graphic
figure(2)
set(gca, 'Color', 'w');
hold on
plot(date(2:end), PIB(2:end), 'k', 'LineWidth', 2)
plot(date(2:end), mean(mu, 2), 'b')
%shade_interval(date(2:end), quantile(mu, 0.16, 2), 'b--', date(2:end), quantile(mu, 0.84, 2), '--b', 'FillType', [1 2; 2 1], 'FillAlpha', 0.1)
legend('Produto','Produto Potencial')
grid on
grid minor

%% f)

% Output gap graphic
figure(3)
set(gca, 'Color', 'w');
hold on
plot(date(2:end), zeros(T-1), 'k-.', 'LineWidth', 1.5)
plot(date(2:end), mean(eta, 2)*100, 'b-', 'LineWidth', 1.5)
%shade_interval(date(2:end), quantile(eta, 0.16, 2), 'b:', date(2:end), quantile(eta, 0.84, 2), 'b:', 'FillType', [1 2; 2 1], 'FillAlpha', 0.1)
grid on
grid minor


