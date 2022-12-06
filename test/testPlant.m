%% HDD driver model
%% tf
k = 3.25e8;
tau = 2.3e-5;
a1 = 1;     zeta1 = 0.6;        omega1 = 30;
a2 = -1.4;  zeta2 = 0.06;       omega2 = 4700;
a3 = 0.4;   zeta3 = 0.02;       omega3 = 6300;
plant1_tf = tf(a1, [1 2*zeta1*omega1 omega1^2]);
plant2_tf = tf(a2, [1 2*zeta2*omega2 omega2^2]);
plant3_tf = tf(a3, [1 2*zeta2*omega3 omega3^2]);
plant_delay_pade_tf = tf([-0.5*tau 1], [0.5*tau 1]);
plant_all_tf = k*plant1_tf*plant2_tf*plant3_tf*plant_delay_pade_tf;
% figure;
% bode(plant_all_tf)
%% ss
plant_all_ss_c = ss(plant_all_tf);
sampling_frequency = 5400;
plant_all_ss_d = c2d(plant_all_ss_c, 1/sampling_frequency);

%% control
% control_num_1 = [1 2900];
% control_num_2 = [1 2544.7];
% control_num_3 = [1 90.6];
% control_num_all = conv(control_num_1, conv(control_num_2, control_num_3));
% control_den_1 = [1 13765.7];
% control_den_2 = [1 2500];
% control_den_3 = [1 0];
% control_den_all = conv(control_den_1, conv(control_den_2, control_den_3));
% control_all_tf = tf(0.17*control_num_all, control_den_all);
% openloop_all_tf = control_all_tf*plant_all_tf;
% figure;
% bode(openloop_all_tf)

%% generate covariance matrix
x_size = 7; y_size = 1; u_size = 1;
covariance_A = zeros(x_size);
covariance_B = zeros(y_size);
covariance_C = zeros(u_size);
while min(eig(covariance_A)) <= 0
    covariance_b1 = rand(x_size);
    covariance_A = covariance_b1*(covariance_b1.');
end
while min(eig(covariance_B)) <= 0
    covariance_b2 = rand(y_size);
    covariance_B = covariance_b2*(covariance_b2.');
end
while min(eig(covariance_C)) <= 0
    covariance_b3 = rand(u_size);
    covariance_C = covariance_b3*(covariance_b3.');
end

covariance = blkdiag(covariance_A, covariance_B, covariance_C);

%% define
load('paraPlant.mat')
para.A = plant_all_ss_d.A;
para.B = plant_all_ss_d.B;
para.C = plant_all_ss_d.C;
para.D = plant_all_ss_d.D;
para.Cov = covariance;

