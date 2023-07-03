function [x_up,P_up,y_f,V_f,x_f,P_f,x0,P0] = kalmanfilter(F,G,Q,H0,H1,R,y,x0,P0)
% -------------------------------------------------------------------------------------
% Command:  kalmanfilter(F,G,Q,H0,H1,R,y,x0,P0)
% Purpose:  Compute the filtered values and one-period ahead from the kalman filter 
% recursions for a constant parameter state-space model as descrived below:
% State-space model with constant parameters
% x(t) = F*x(t-1) + G*v(t)  , v(t)~N(0,Q) iid
% y(t) = H0 + H1*x(t) + w(t), w(t)~N(0,R) iid
% Inputs:   F  -- r x r matrix
%           G  -- r x s matrix with s<=r(the matrix G allows less shocks than states)
%           Q  -- s x s covariance matrix
%           H0 -- m x 1 vector of constants
%           H1 -- m x s matrix
%           R  -- m x m covariance matrix
%           y  -- T x m data matrix
%           x0 (optional) -- r x 1 initial state vector: x(t|t)
%           P0 (optional) -- r x r initial state covariance matrix: x(t|t)
% Outputs:  x_up -- T x r filtered values for states
%           P_up -- r x r x T filtered covariance matrix for states 
%           y_f (optional) -- T x m one-period ahead forecast for observables: y(t|t-1)for t=1,...,T
%           V_f (optional) -- m x m x T one-period ahead MSE for observables: V(t|t-1) for t=1,...,T
%           x_f (optional) -- T x r one-period ahead forecast for states: x(t+1|t) for t=1,...,T
%           P_f (optional) -- r x r x T one-period ahead MSE for states: P(t+1|t) for t=1,...,T
%           x0  (optional) -- r x 1 initial state vector
%           P0  (optional) -- r x r initial state covariance matrix
%           for t=1,...,T where T is the size of y
% -------------------------------------------------------------------------------------
% Date: 23/06/2020
% Author: Marcel Ribero
% If you find any errors, please contact marcel.ribeiro@fgv.br
% -------------------------------------------------------------------------------------

% Initialize the State Vector at the Stationary Distribution
y = transpose(y);
[n,T]    = size(y);
r        = size(F,1);
x_up     = zeros(r,T);   %x(t|t)
P_up     = zeros(r,r,T); %P(t|t)
	if nargout > 2 %If asked for the filtered values only, does not need to compute all objects
        x_f      = zeros(T,r);   %x(t+1|t) for t=1,...,T
        P_f      = zeros(r,r,T); %P(t+1|t) for t=1,...,T
        y_f      = nan(T,n);     %y(t|t-1) for t=1,...,T 
        V_f      = zeros(n,n,T); %V(t|t-1) for t=1,...,T
	end
    
    QQ = G*Q*transpose(G); % QQ is the state covariance matrix of the compound error e(t) = G*v(t)
    if nargin < 8  % Use unconditional mean if not provided
        x0 = zeros(r,1);     %Initial value as unconditional mean for x(0|0)
        % Initial values for of the covariance matrix of the state P(0|0)
        P0        = inv(eye(r*r) - kron(F,F))*reshape(QQ,[],1);
        P0 = reshape(P0,r,r);
        % Initialization is faster with dlyap (if toolbox is available)
        % P_up(:,:,1) = dlyap(F,Q_all);
    end
    % Initial values for x(1|0) and P(1|0)
    xf = F*x0;
    Pf = F*P0*transpose(F) + QQ;
    
%% Kalman Filter Recursions
for t = 1:T 
    %% Measurement forecast step (y(t+1|t) and V(t+1|t))
    yf = H0 + H1*xf;
    Vf = H1*Pf*transpose(H1) + R;
                
    %% Update step 
    Kgain = Pf*transpose(H1)/Vf;
    x_up(:,t) = xf + Kgain*(y(:,t) - yf);
    P_up(:,:,t) = Pf - Kgain*H1*Pf;

	%% State forecast step (x(t+1|t) and P(t+1|t))
    xf = F*x_up(:,t);
    Pf = F*P_up(:,:,t)*transpose(F) + QQ;

    if nargout > 2 % if requesting more than filtered values
        y_f(t,:)   = yf';  %y_f refers to y(t|t-1) where 't' refers to the loop index 
        P_f(:,:,t) = Pf;   %P_f refers to P(t+1|t) where 't' refers to the loop index
        x_f(t,:)   = xf';  %x_f refers to x(t+1|t) where 't' refers to the loop index
        V_f(:,:,t) = Vf;   %V_f refers to V(t|t-1) where 't' refers to the loop index
    end
end
    x_up = transpose(x_up);

end
