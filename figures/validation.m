%{
Run the validation testcases
%}

clear all; close all; clc;

% run_valid('2d_UU_UU-type=0-orderdiff=2',2);
% run_valid('2d_UU_UU-type=0-orderdiff=3',4);

% run_valid('2d_UU_UE-type=0-orderdiff=2',2);
% run_valid('2d_UU_UE-type=0-orderdiff=3',4);

% run_valid('2d_UU_UO-type=0-orderdiff=2',2);
run_valid('2d_UU_UO-type=0-orderdiff=3',4);


run_valid('2d_UU_OU-type=0-orderdiff=3',4);

