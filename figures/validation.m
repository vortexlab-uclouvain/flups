%{
Run the validation testcases
%}

clear all; close all; clc;


%------ unb unb - unb unb - unb unb
% run_valid('3d_444444_typeGreen=0',2); 
% run_valid('3d_444444_typeGreen=2',2);
% run_valid('3d_444444_typeGreen=3',4);
% run_valid('3d_444444_typeGreen=4',6);

%------ even unb - unb unb - unb unb
% run_valid('3d_044444_typeGreen=0',2); 
% run_valid('3d_044444_typeGreen=2',2);
% run_valid('3d_044444_typeGreen=3',4);
% run_valid('3d_044444_typeGreen=4',6);

%------ unb unb - unb odd - unb unb
% run_valid('3d_444144_typeGreen=0',2); 
% run_valid('3d_444144_typeGreen=2',2);
% run_valid('3d_444144_typeGreen=3',4);
% run_valid('3d_444144_typeGreen=4',6);

%------ unb unb - unb unb - unb odd
% run_valid('3d_444441_typeGreen=0',2); 
run_valid('3d_444441_typeGreen=2',2);
% run_valid('3d_444441_typeGreen=3',4);
% run_valid('3d_444441_typeGreen=4',6);

% run_valid('3d_UU_UU_UU-type=0-typeGreen=2',2);
% run_valid('3d_UU_UU_UU-type=0-typeGreen=3',4);

% run_valid('3d_EU_UU_UU-type=0-typeGreen=2',2);
% run_valid('3d_OU_UU_UU-type=0-typeGreen=2',2);
% run_valid('3d_UU_OU_UU-type=0-typeGreen=2',2);

% run_valid('2d_UU_UU-type=0-orderdiff=2',2);
% run_valid('2d_UU_UU-type=0-orderdiff=3',4);
% 
% run_valid('2d_UU_UE-type=0-orderdiff=2',2);
% run_valid('2d_UU_UE-type=0-orderdiff=3',4);
% run_valid('2d_UU_EU-type=0-orderdiff=2',2);
% run_valid('2d_UU_EU-type=0-orderdiff=3',4);
% run_valid('2d_UE_UU-type=0-orderdiff=2',2);
% run_valid('2d_UE_UU-type=0-orderdiff=3',4);
% run_valid('2d_EU_UU-type=0-orderdiff=2',2);
% run_valid('2d_EU_UU-type=0-orderdiff=3',4);
% 
% 
% run_valid('2d_UU_UO-type=0-orderdiff=2',2);
% run_valid('2d_UU_UO-type=0-orderdiff=3',4);
% run_valid('2d_UU_OU-type=0-orderdiff=2',2);
% run_valid('2d_UU_OU-type=0-orderdiff=3',4);
% run_valid('2d_OU_UU-type=0-orderdiff=2',2);
% run_valid('2d_OU_UU-type=0-orderdiff=3',4);
% run_valid('2d_UO_UU-type=0-orderdiff=2',2);
% run_valid('2d_UO_UU-type=0-orderdiff=3',4);
