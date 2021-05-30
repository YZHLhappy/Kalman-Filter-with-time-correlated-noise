function [Pp,Pv,Kp,Kv,p,v] = KF_Perturbed_P(n,l_obs,A,x_ini,P_ini,phi,psi,R,Q,B_dis)
% Perturbed-P kalman filtering: from Kedong Wang, Yong Li,Chris Rizos article (2012)
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

x_plus = x_ini;
H = A;
PHIa = [phi,[0;0];0,0,psi];

GAMMAa = [B_dis,[0;0];0,0,1];
Ha = [H,1];
Qa = GAMMAa * [Q,[0;0];0,0,R]*GAMMAa';
Pa_plus = [P_ini,[0;0];0,0,sigma_ce^2];

lambada = 10^-6;                                                           % lambada is a small value (such as 10^-6)
xa_plus = [x_plus;0];

for i = 1:n-1   
    xa_minus = PHIa * xa_plus ;   
    Pa_minus = PHIa * Pa_plus * PHIa' + Qa;
    Ka = Pa_minus * Ha'* inv(Ha * Pa_minus * Ha');
    xa_plus = xa_minus + Ka * (l_obs(i) - Ha * xa_minus);
    Pa_plus = Pa_minus - Ka * Ha * Pa_minus + lambada * eye(3);

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





















