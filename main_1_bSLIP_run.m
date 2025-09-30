% Inverted Pendulum Walking Simulation

clear;
close all;
clc;

% define parameters, initial conditions, time settings
setup_inverted_pendulum_walk;

out = sim("inverted_pendulum_walk.slx");
