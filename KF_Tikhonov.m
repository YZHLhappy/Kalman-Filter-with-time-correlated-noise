function [Pp,Pv,Kp,Kv,p,v] = KF_Tikhonov(n,l_obs,A,x_ini,P_ini,phi,psi,R,Q,B_dis)
% Tikhonv kalman filtering: from Kedong Wang, Yong Li,Chris Rizos article (2012)
% "Practical Approaches to Kalman Filtering with Time-Correlated Measurement Errors"
% Add an a after the parameter to indicate augmentation

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
% B_dis: the transition matrix of the noise part

% Output: 
% Pp: P-Matrix(1,1),position-term 
% Pv: P-Matrix(2,2),velocity-term
% Kp: K-Vector(1,1),position-term
% Kv: K-Vector(1,1),velocity-term
% p:position; v:velocity

H = A;
PHIa = [phi,[0;0];0,0,psi];                     
GAMMAa = [B_dis,[0;0];0,0,1];
Ha = [H,1];
Qa = GAMMAa*[Q,[0;0];0,0,R]*GAMMAa';
Pa_plus_standard = [P_ini,[0;0];0,0,sigma_ce^2];

N = 50;                                                                    % N is a threshold iteration number
Beta = 1*10^-4;                                                            % Beta: a small positive value,such as:1*10^-4;
[m,~] = size(Ha);
y = ones(m,1);
epsilon = 1*10^-6;                                                         % is a small value,such as:1*10^-6
xa_plus = [x_ini;0];

for i = 1:n-1
    xa_minus = PHIa * xa_plus;
    Pa_minus = PHIa * Pa_plus_standard * PHIa' + Qa;
    
    W_k = Ha * Pa_minus * Ha';
    u_alpha_n = inv(W_k'*W_k)*W_k'*y;
    alpha_n = Beta * (1+u_alpha_n' * u_alpha_n)-(y-W_k*u_alpha_n)'*(y-W_k*u_alpha_n)/(1+u_alpha_n' * u_alpha_n);
    for j = 1:N
        u_alpha_nplus1 = inv(W_k'*W_k + alpha_n * eye(1))*W_k'*y;
        alpha_nplus1 = Beta * (1+u_alpha_nplus1' * u_alpha_nplus1)-(y-W_k*u_alpha_nplus1)'*(y-W_k*u_alpha_nplus1)/(1+u_alpha_nplus1' * u_alpha_nplus1);
        norm_u_alpha = norm(u_alpha_nplus1 - u_alpha_n);
        if norm_u_alpha < epsilon
            break
        else 
            u_alpha_n = u_alpha_nplus1;
            alpha_n = alpha_nplus1;
        end
    end
    alpha = alpha_nplus1;                                                  % regularization parameter
    
    Ka_standard = Pa_minus* Ha'*inv(alpha * eye(1) + W_k'*W_k)*W_k';       % kalman gain
    xa_plus = xa_minus + Ka_standard * (l_obs(i) - Ha * xa_minus);
    Pa_plus = Pa_minus - Ka_standard * Ha * Pa_minus;
    
    sigma_p(i) = sqrt(Pa_plus(1,1));                                       % standard deviation position
    sigma_v(i) = sqrt(Pa_plus(2,2));                                       % standard deviation velocity       
    Kp(i) = Ka_standard(1,1);                                              % output; kalman_position   
    Kv(i) = Ka_standard(2,1);                                              % output; kalman_position
    p(i)= xa_plus(1);                                                      % output: position
    v(i)= xa_plus(2);                                                      % output: velocity
end
Pp = sigma_p';                                                             % output: standard deviation_position
Pv = sigma_v';                                                             % output: standard deviation_velocity
end














