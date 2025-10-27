% gravity & masses/lengths
g = 9.81;
m_foot = 1;
r_foot = 0.01;
m_total = 80;
m = m_total - 2*m_foot;

l_total = 1.0;
l_nominal = l_total - r_foot;
l_TD = l_nominal;

% ground
plane_height = 0.1; 
plane_x = 10;
plane_z = 1;
fric_mu_s = 1.0;
fric_mu_d = 1.0;
contact_k = 1e6;
contact_c = 1e6;

% initial COM state 
x0 = 0;  
y0 = 1;  
dx0 = 4; 
dy0 = 0;

% initial legs
phi0_L = 20*pi/180;
phi0_R = -phi0_L;
dphi0_L = 0;
dphi0_R = 0;
l0_L = l_nominal;
l0_R = l_nominal;
dl0_L = 0; 
dl0_R = 0;

% stance passive spring-damper
k_leg = 20000;
c_leg = 50;

% swing controls (hip PD + leg-length PD)
phi_TD_deg = 20;
phi_TD = phi_TD_deg*pi/180;
hip_P = 200; 
hip_D = 10;
leg_P = 5000;
leg_D = 100;

% "short leg" in swing to avoid toe drag
l_swing_short = 0.5;
phi_th = 1*pi/180;

% termination conditions
body_ground_clearance_stop = 0.0;
max_double_support_duration = 1.0;
max_leg_collapse_in_stance = 0.5;
max_ground_penetration = 0.05;

% Raibert touchdown params
phi_bias = deg2rad(20);
dx_ref = 5.0;
k_v = 0.11;
phi_min = deg2rad(10); 
phi_max = deg2rad(28); 


% visuals & sim time
l1 = 0.4; l2 = 0.4;
t_end = 5;