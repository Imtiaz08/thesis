clear all; clc; close all;

addpath external_tools
addpath external_tools/GPS_Broadcast_Orbits
addpath external_tools/danipascual-GNSS-matlab-f2d8287
addpath external_tools/danipascual-GNSS-matlab-f2d8287/prn_codes
addpath external_tools/GPScode
addpath external_tools/genE5code
addpath external_tools/EllipsoidPlaneIntersection
addpath external_tools/hline_vline

addpath functions



set(0,'defaultaxesfontsize', 16);
set(0,'defaultlinelinewidth', 2);

dataFolder = '../Data/'; % General data folder

featureFolder = [dataFolder 'Features/'];
figureFolder = [dataFolder 'Figures/'];
labsatFolder = [dataFolder 'Labsat/'];
matFolder = [dataFolder 'Mat/'];
NMEAfolder = [dataFolder 'NMEA/'];
tripFolder = [dataFolder 'Triplog/'] ;

% Define constants
const.c = 299792458; % Speed of light, m/s
const.eps0 = 8.854187817e-12; % Permeability in vacuum, F/m
const.R_e = 6371008.7714; %6.3781E6; % Earth radius, m

const.LS3.f_samp = 16.368e6; % Hz

const.GPS.L1.CA.f_carrier = 1.57542e9; % Hz
const.GPS.L1.CA.T_code = 1e-3; % s
const.GPS.L1.CA.L_code = 1023; % bits
const.GPS.L1.CA.f_chip = 1.023e6; % bits/s
const.GPS.L1.CA.T_data = 50e-3; % s

const.GPS.L2.CL.f_carrier = 1.2276e9; 
const.GPS.L2.CL.T_code = 1.5;
const.GPS.L2.CL.L_code = 767250; 
const.GPS.L2.CL.f_chip = 0.5115e6; 

const.GPS.L2.CM.f_carrier = 1.2276e9; 
const.GPS.L2.CM.T_code = 20e-3;
const.GPS.L2.CM.L_code = 10230; 
const.GPS.L2.CM.f_chip = 0.5115e6; 

const.GPS.L5.SOL.f_carrier = 1176.45e6; % Secondary codes?????
const.GPS.L5.SOL.T_code = 1e-3;
const.GPS.L5.SOL.L_code = 10230; 
const.GPS.L5.SOL.f_chip = 10.23e6; 

const.GAL.E1.OS.f_carrier = 1.57542e9;
const.GAL.E1.OS.T_code = 4e-3;
const.GAL.E1.OS.L_code = 4092;
const.GAL.E1.OS.f_chip = 1.023e6;
const.GAL.E1.OS.f_subcarrier = 1.023e6; 
const.GAL.E1.OS.m = 6;
const.GAL.E1.OS.n = 1;
const.GAL.E1.OS.alpha = sqrt(10/11);
const.GAL.E1.OS.beta = sqrt(1/11);

const.GAL.E5.OS.f_carrier = 1191.795e6;
const.GAL.E5.OS.T_code = 1e-3;
const.GAL.E5.OS.L_code = 10230;
const.GAL.E5.OS.f_chip = 15.345e6;
const.GAL.E5.OS.f_subcarrier = 15.345e6;

const.GAL.E5a.OS.f_carrier = 1176.45e6; 
const.GAL.E5a.OS.T_code = 1e-3;
const.GAL.E5a.OS.L_code = 10230;
const.GAL.E5a.OS.f_chip = 15.345e6;
const.GAL.E5a.OS.f_subcarrier = 15.345e6;

const.GAL.E5b.OS.f_carrier = 1207.14e6;
const.GAL.E5b.OS.T_code = 1e-3;
const.GAL.E5b.OS.L_code = 10230;
const.GAL.E5b.OS.f_chip = 15.345e6;
const.GAL.E5b.OS.f_subcarrier = 15.345e6;

%Source: Hann01, p. 37
% Material constants are in reality also frequency-dependent!
const.cond.concrete = 2e-5; % Sigma, conductivity
const.cond.dryGround = 1e-5;
const.cond.medDryGround = 4e-2;
const.cond.wetGround = 2e-1;
const.cond.freshWater = 2e-1;
const.cond.seaWater = 4;
const.cond.glass = 0.0025;
const.cond.steel = 1.45e6;

const.relPerm.concrete = 3; % Epsilon, relative permittivity
const.relPerm.dryGround = 4;
const.relPerm.medDryGround = 7;
const.relPerm.wetGround = 30;
const.relPerm.freshWater = 80;
const.relPerm.seaWater = 20;
const.relPerm.glass = 7;
const.relPerm.steel = 1;
