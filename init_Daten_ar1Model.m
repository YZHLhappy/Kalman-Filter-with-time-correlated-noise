function [n,tau,Ts,n_simulation,P_ini,x_ini,velocity,A,phi,std_Q_v,...
    std_Q_p,psi,sigma_ce,Q] = init_Daten_ar1Model(situation_Q)
% This function contains initial parameters.
% Input:
% situation_Q = 0: without process noise Q
% situation_Q = 1: with process noise Q

% Output:
% n,tau,Ts,n_simulation,P_ini,x_ini,velocity,A,phi,std_Q_v,std_Q_p,psi,sigma_ce,Q
n = 1000;                                                                  % length of the observation vector
tau = 10;                                                                  % time constant of correlated errors
Ts = 1;                                                                    % time increment 
psi = exp(-Ts/tau);                                                        % time-correlated coefficient
sigma_ce = 5;                                                              % Standard deviation of correlated errors,1/5/10...
n_simulation = 100;                                                        % monte carlo simulation times

% MODEL: observation: l_obs=[velocity], state: x=[position;velocity]
P_ini = [10 0; 0 10];                                                      % initial value of P
x_ini = [0;1];                                                             % initial value of state
velocity = 1;                                                              % 1[ms-1]
A = [0 1];                                                                 % measurement matrix
phi_constinous = [0 1;0 0];                                                % design matrix continous case:constant velocity
phi = expm(phi_constinous);                                                % discrete design matrix

% Q: process noise matrix
if situation_Q == 0                                                        % without process noise Q
    factQ = 0;
else
    factQ = 1;                                                             % with process noise Q
end  
std_Q_v = 1e-1;                                                            % standard deviation of Qv
std_Q_p = std_Q_v^2;                                                       % standard deviation of Qp
B_Q = [std_Q_p^2,0;0,std_Q_v^2];
Am = [-phi_constinous*Ts B_Q*B_Q'*Ts;zeros(2) phi_constinous'*Ts];
Br_2 = expm(Am);
B_disQ_2 = phi*Br_2(1:2,3:4);
Q = factQ*B_disQ_2;
end

