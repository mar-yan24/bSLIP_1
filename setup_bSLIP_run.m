% model parameters
g = 9.81; % []
m_foot = 1;
r_foot = 0.01;
m = 80 - 2*m_foot; % [kg]
l_total = 1;
l = l_total-r_foot; % [m]

% ground
plane_height = 0.1;
plane_x = 10;
plane_z = 1;

% set main initial conditions
x0 = 1;
v = 1.3;
phi_max = 15*pi/180; % rad

% dervie leg initial condition
phi0_L = phi_max;
d_phi0_L = -v/l_total; % rad/sec
phi0_R = -phi0_L; % rad
d_phi0_R = -d_phi0_L; % rad/sec

% dervie COM initial condition
d_x0 = v*cos(phi0_L); % [m/s]
y0 = l_total*cos(phi0_L);
d_y0 = v*sin(phi0_L);

% swing leg control
l_min = 0.8;
hip_P = 200; % P gain for hip angle control
hip_D = 10; % D gain for hip angle control

phi_th = 1*pi/180;

% visualization
l1 = 0.4; % length of top leg segment
l2 = 0.4; % length of lower leg segment

% % Time settings
t_end = 3;
