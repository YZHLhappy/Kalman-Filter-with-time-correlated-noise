n = 1000;
psi = 0.4;
sigma_ce = 5;
Ts = 1;
tau = 10;
v = 1;
[u_ar1,l_obs,R] = Observation_with_AR1Noise_withoutQ(n,psi,v,sigma_ce,Ts,tau);

subplot(4,1,1)
plot(u_ar1)
ylim([-30,30])
title('AR(1) Rauschen,sigma(white noise)=5,phi1=0.4')

phi_AR2 = [0.4,0.2];
sigma_w = sigma_ce;
[u_ar2,l_obs_ar2,R] = Observation_with_AR2Noise_withoutQ(n,phi_AR2,v,sigma_w);
subplot(4,1,2)
ylim([-30,30])
plot(u_ar2)
title('AR(2) Rauschen,sigma(white noise)=5,phi1=0.4,phi2=0.2')

phi_AR4 = [0.4,0.2,0.1,0.05];
phi = [1,1;0,1];    
H = [0,1];
[u_ar4,l_obs_ar4,R] = Observation_with_AR4Noise_withoutQ(n,phi_AR4,v,sigma_w,phi,H);
subplot(4,1,3)
plot(u_ar4)
ylim([-30,30])
title('AR(4) Rauschen,sigma(white noise)=5,phi1=0.4,phi2=0.2,phi3=0.1,phi4=0.05')

A=0.8;
[u_out,l_obs,R] = ffgn_Noise_Observation(n,v,A);
subplot(4,1,4)
plot(u_out)
ylim([-30,30])
title('Fractional Gaussian noise,sigma=5,H=0.8,force=0')














