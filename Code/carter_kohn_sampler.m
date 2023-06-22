function xdraws = carter_kohn_sampler(F,G,Q,H0,H1,R,y,x0,P0)
% ---------------------------------------------------------------------------------
% Command:  carter_kohn_sampler(F,G,Q,H0,H1,R,y,x0,P0)
% Purpose:  Ramdomly draw the series of states from a constant parameter state-space
% model as descrived below:
% State-space model with constant parameters
% x(t) = F*x(t-1) + G*v(t)  , v(t)~N(0,Q) iid
% y(t) = H0 + H1*x(t) + w(t), w(t)~N(0,R) iid
% Inputs:   F  -- r x r matrix
%           G  -- r x s matrix
%           Q  -- s x s covariance matrix
%           H0 -- m x 1 vector of constants
%           H1 -- m x s matrix
%           R  -- m x m covariance matrix
%           y  -- T x m data matrix
%           x0 -- r x 1 initial state vector
%           P0 -- r x r initial state covariance matrix
% Outputs:  xdraws T+1 x r random drwas
% ----------------------------------------------------------------------------------
% Date: June/2020
% Author: Marcel Ribero
% If you find any errors, please contact marcel.ribeiro@fgv.br
% ----------------------------------------------------------------------------------

% Run Kalman filter to get x(t|t) and P(t|t)
[x_up,P_up,~,~,x_f,P_f] = kalmanfilter(F,G,Q,H0,H1,R,y,x0,P0);

% Define objects
[T,n_state] = size(x_up);
xdraws = zeros(T,n_state);

% Draw x_T from N(x(T|T),P(T|T))
xdraws(T,:) = mvnrnd(x_up(T,:),P_up(:,:,T),1);

% Calculate backwards the smoothed values and smoothed MSE
for t=T-1:-1:1
    J= P_up(:,:,t)*transpose(F)/P_f(:,:,t);           % Compute the "Smoother gain"
    mu = x_up(t,:)' + J*(xdraws(t+1,:)' - x_f(t,:)'); % mean of x(t)|Y(t),x_(t+1)
    cov = P_up(:,:,t) - J*F*P_up(:,:,t);              % cov  of x(t)|Y(t),x_(t+1)
    xdraws(t,:) = mvnrnd(mu,cov,1);                   % Draw x_t from N(mean,cov)
end

end








