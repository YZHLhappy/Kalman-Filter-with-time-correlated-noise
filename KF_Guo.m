function [Pp,Pv,Kp,Kv,p,v] = KF_Guo(n,l_obs,A,x_ini,P_ini,phi,psi,R,Q)
% Guo kalman filtering: from Guobin Chang's article (2014)
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

x_plus = x_ini;
P_plus = P_ini;
H = A;

z_k = 0;                                                                   % The difference between the observations before and after the moment
for j = 1:LE_simu-1
    z_k(j+1) = l_obs(j+1) - psi*l_obs(j);
end

for i = 1:n-1  
    x_minus = phi * x_plus;
    P_minus = phi * P_plus * phi' + Q;
    P_ = phi * P_plus;
    n_Guo = z_k(i) - H * x_minus + psi_mat * H * x_plus;
    Sigma_Guo = H * P_minus * H'+psi_mat * H * P_plus * H'* psi_mat'+ R...
        -H * phi * P_plus * H'* psi_mat'-psi_mat * H * P_plus * phi'* H';
    P_hat = P_minus * H'- P_ * H'* psi_mat';
    K_Guo = P_hat * inv(Sigma_Guo);                                        % kalman gain
    x_plus = x_minus + K_Guo * n_Guo;
    P_plus = P_minus - K_Guo * Sigma_Guo * K_Guo';

    sigma_p(i) = sqrt(P_plus(1,1));                                        % standard deviation: position
    sigma_v(i) = sqrt(P_plus(2,2));                                        % standard deviation: velocity  
    Kp(i) = K_Guo(1,1);                                                    % output: kalman gain_position
    Kv(i) = K_Guo(2,1);                                                    % output: kalman gain_velocity
    p(i)=x_plus(1);                                                        % output: position
    v(i)=x_plus(2);                                                        % output: velocity
end
Pp = sigma_p';                                                             % output: standard deviation_position
Pv = sigma_v';                                                             % output: standard deviation_velocity
end
















