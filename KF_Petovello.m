function [Pp,Pv,Kp,Kv,p,v] = KF_Petovello(n,l_obs,A,x_ini,P_ini,phi,psi,R,Q)
% Petovello kalman filtering from Guobin Chang's article (2014)
% "On kalman filter for linear system with colored measurement noise"

% Input:
% n: the length of the simulation time
% l_obs: observations
% A: observation matrix
% x_ini: initial state vector
% P_ini: initial error matrix
% phi: state transition matrix
% psi: time-correlated noise coefficient
% R: observation noise matrix
% Q：system noise matrix

% Output: 
% Pp: P-Matrix(1,1),position-term 
% Pv: P-Matrix(2,2),velocity-term
% Kp: K-Vector(1,1),position-term
% Kv: K-Vector(1,1),velocity-term
% p:position; v:velocity

H = A;
Gamma = psi * H * inv(phi);
L = H - Gamma;
R_ = Gamma * Q * Gamma'+ R;
S_ = Q * Gamma';
P_plus = P_ini;
x_plus = x_ini;

z_k = 0;                                                                   % The difference between the observations before and after the moment
for j = 1:LE_simu-1
    z_k(j+1) = l_obs(j+1) - psi_mat.*l_obs(j);
end

for i = 1:n-1  
    x_minus = phi * x_plus;
    P_minus = phi * P_plus * phi' + Q;
    d = z_k(i) - L * x_minus;
    sigma = L * P_minus * L'+ R_ - L * S_ - S_'* L';
    P_xhat = P_minus * L'+ S_; 
    K = P_xhat * inv(sigma);                                               % kalman gain
    x_plus = x_minus + K * d;
    P_plus = P_minus - K * sigma * K';

    sigma_p(i) = sqrt(P_plus(1,1));                                        % standard deviation: position
    sigma_v(i) = sqrt(P_plus(2,2));                                        % standard deviation: velocity  
    Kp(i) = K(1,1);                                                        % output: kalman gain_position
    Kv(i) = K(2,1);                                                        % output: kalman gain_velocity
    p(i)=x_plus(1);                                                        % output: position
    v(i)=x_plus(2);                                                        % output: velocity
end
Pp = sigma_p';                                                             % output:standard deviation_position
Pv = sigma_v';                                                             % output:standard deviation_velocity
end



