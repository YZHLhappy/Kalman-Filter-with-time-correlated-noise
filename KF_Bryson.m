function [Pp,Pv,Kp,Kv,p,v] = KF_Bryson(n,l_obs,A,x_ini,P_ini,phi,psi,R,Q)
% Brsyon kalman filtering from Guobin Chang's article (2014)
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
x_plus = x_ini;
P_plus = P_ini;

z_k = 0;                                                                   % The difference between the observations before and after the moment
for j = 1:n-1
    z_k(j+1) = l_obs(j+1) - psi*l_obs(j);
end
G = H * phi - psi_mat * H;
R_ = H * Q * H' + R;  
sigma = G * P_plus * G' + R_;
S = Q * H';
J = S * inv(R_);
F = phi - J * G;
Q_plus = Q - J * R_ * J';
for i = 1:LE_simu-1
    P_minus = P_plus * G'; 
    K = P_minus * inv(sigma);                                              % kalman gain
    x_plus = (F - F * K * G) * x_plus + (F * K + J) * z_k(i);
    P_plus = F * P_plus * F'- F * K * sigma * K' * F' + Q_plus; 

    sigma_p(i) = sqrt(P_plus(1,1));                                        % standard deviation: position 
    sigma_v(i) = sqrt(P_plus(2,2));                                        % standard deviation: velocity
    p_est(i)=x_plus(1);                                                    % position 
    v_est(i)=x_plus(2);                                                    % velocity
    Kp(i) = K(1,1);                                                        % kalman gain: position
    Kv(i) = K(2,1);                                                        % kalman gain: velocity 
end
Pp = sigma_p';                                                             % output:standard deviation: position
Pv = sigma_v';                                                             % output:standard deviation: velocity
p = p_est';                                                                % output: position
v = v_est';                                                                % output: velocity
end








