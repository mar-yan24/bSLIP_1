% Gravity & masses/lengths
g = 9.81;
m_foot = 1;
r_foot = 0.01; % foot radius
m_total = 80; % body+feet
m = m_total - 2*m_foot;  % point-mass body used by Simscape "body" block

l_total = 1.0; % hip-to-foot when straight down
l_nominal = l_total - r_foot; % "rest" (effective) leg length used in stance spring (a.k.a. l_min in table)
l_TD = l_nominal; % target leg length at touchdown (leading swing reaches this before contact)

% Ground (world frame)
plane_height = 0.1; 
plane_x = 10;
plane_z = 1;
fric_mu_s = 1.0; % static friction
fric_mu_d = 1.0; % dynamic friction
contact_k = 1e6;
contact_c = 1e6;

% Initial COM state (flight)
x0 = 0;  
y0 = 1;  
dx0 = 4; 
dy0 = 0;

% Initial legs
phi0_L = 20*pi/180;   % hip angle: 0 = down, + = flexion
phi0_R = -phi0_L;
dphi0_L = 0;
dphi0_R = 0;
l0_L = l_nominal;
l0_R = l_nominal;
dl0_L = 0; 
dl0_R = 0;

% Stance passive spring-damper (bSLIP)
k_leg = 20000;          % N/m
c_leg = 50;             % N*s/m

% Swing controls (hip PD + leg-length PD)
phi_TD_deg = 20; % touchdown hip angle (used as "leading-swing" target)
phi_TD = phi_TD_deg*pi/180;
hip_P = 200; % PD gains for hip angle control
hip_D = 10;
leg_P = 5000; % PD gains for leg length control in swing
leg_D = 100;

% "short leg" in swing to avoid toe drag (trailing & single swing)
l_swing_short = 0.5; % fast PD toward this when the leg is trailing/single swing
phi_th = 1*pi/180; % small angle threshold

% Termination conditions (check in Simulink using Assert/Stop blocks)
body_ground_clearance_stop = 0.0;    % stop if body y <= 0
max_double_support_duration = 1.0;    % s
max_leg_collapse_in_stance = 0.5;    % m (i.e., l_nominal - l > 0.5 => stop)
max_ground_penetration = 0.05;   % m

% Visuals & sim time
l1 = 0.4; l2 = 0.4;
t_end = 5;