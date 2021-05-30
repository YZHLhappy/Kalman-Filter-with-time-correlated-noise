function [Pp,Pv,Kp,Kv,p,v] = KF_AR2(n,l_obs,A,x_ini,P_ini,phi,R,Q,psi_AR2)
% state augmentation approach：for AR(2)-noise
% Add an a after the parameter to indicate augmentation

% Input:
% n: the length of the simulation time
% l_obs: observations
% A: observation matrix
% x_ini: initial state vector
% P_ini: initial error matrix
% phi: state transition matrix
% psi_AR2: time-correlated noise coefficient
% R: observation noise matrix
% Q：system noise matrix

% Output: 
% Pp: P-Matrix(1,1),position-term 
% Pv: P-Matrix(2,2),velocity-term
% Kp: K-Vector(1,1),position-term
% Kv: K-Vector(1,1),velocity-term
% p:position; v:velocity

GAMMA =[1,0;0,1];                                                          %The noise transfer matrix here is a fixed value. Need to be modified according to the model.
GAMMAa = [GAMMA,[0,0;0,0];[0,0;0,0],[1,0;0,1]];
Qa = GAMMAa * [Q,[0,0;0,0];[0,0;0,0],[R,0;0,0]] * GAMMAa';
H = A;
Ha = [H,1,0];
phi_a = [phi,[0,0;0,0];[0,0;0,0],[phi_AR2(1),phi_AR2(2);1,0]];
Pa_plus = [P_ini,[0,0;0,0];[0,0;0,0],[1,0;0,1]];
xa_plus = [x_ini;0;0];

for i = 1:n
    xa_minus = phi_a * xa_plus;
    Pa_minus = phi_a * Pa_plus * phi_a' + Qa;
    Ka_standard = Pa_minus* Ha'*inv(Ha*Pa_minus*Ha');
    xa_plus = xa_minus + Ka_standard * (l_obs(i) - Ha * xa_minus);
    Pa_plus = (eye(4)-Ka_standard * Ha) * Pa_minus;

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

