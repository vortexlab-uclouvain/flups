%{
Run the validation testcases
%}

clear all; close all; clc;


%------ unb unb - unb unb - unb unb
run_valid('3d_444444_typeGreen=0',2); 
run_valid('3d_444444_typeGreen=2',2);
run_valid('3d_444444_typeGreen=3',4);
run_valid('3d_444444_typeGreen=4',6);

%------ x unb - unb unb - unb unb
run_valid('3d_044444_typeGreen=0',2); 
run_valid('3d_044444_typeGreen=2',2);
run_valid('3d_044444_typeGreen=3',4);
run_valid('3d_044444_typeGreen=4',6);

run_valid('3d_144444_typeGreen=0',2); 
run_valid('3d_144444_typeGreen=2',2);
run_valid('3d_144444_typeGreen=3',4);
run_valid('3d_144444_typeGreen=4',6);

%------ unb x - unb unb - unb unb
run_valid('3d_404444_typeGreen=0',2); 
run_valid('3d_404444_typeGreen=2',2);
run_valid('3d_404444_typeGreen=3',4);
run_valid('3d_404444_typeGreen=4',6);

run_valid('3d_414444_typeGreen=0',2); 
run_valid('3d_414444_typeGreen=2',2);
run_valid('3d_414444_typeGreen=3',4);
run_valid('3d_414444_typeGreen=4',6);

%------ unb unb - x unb - unb unb
run_valid('3d_440444_typeGreen=0',2); 
run_valid('3d_440444_typeGreen=2',2);
run_valid('3d_440444_typeGreen=3',4);
run_valid('3d_440444_typeGreen=4',6);

run_valid('3d_441444_typeGreen=0',2); 
run_valid('3d_441444_typeGreen=2',2);
run_valid('3d_441444_typeGreen=3',4);
run_valid('3d_441444_typeGreen=4',6);

%------ unb unb - unb x - unb unb
run_valid('3d_444044_typeGreen=0',2); 
run_valid('3d_444044_typeGreen=2',2);
run_valid('3d_444044_typeGreen=3',4);
run_valid('3d_444044_typeGreen=4',6);

run_valid('3d_444144_typeGreen=0',2); 
run_valid('3d_444144_typeGreen=2',2);
run_valid('3d_444144_typeGreen=3',4);
run_valid('3d_444144_typeGreen=4',6);

%------ unb unb - unb unb - x unb
run_valid('3d_444404_typeGreen=0',2); 
run_valid('3d_444404_typeGreen=2',2);
run_valid('3d_444404_typeGreen=3',4);
run_valid('3d_444404_typeGreen=4',6);

run_valid('3d_444414_typeGreen=0',2); 
run_valid('3d_444414_typeGreen=2',2);
run_valid('3d_444414_typeGreen=3',4);
run_valid('3d_444414_typeGreen=4',6);

%------ unb unb - unb unb - unb x
run_valid('3d_444440_typeGreen=0',2); 
run_valid('3d_444440_typeGreen=2',2);
run_valid('3d_444440_typeGreen=3',4);
run_valid('3d_444440_typeGreen=4',6);

run_valid('3d_444441_typeGreen=0',2); 
run_valid('3d_444441_typeGreen=2',2);
run_valid('3d_444441_typeGreen=3',4);
run_valid('3d_444441_typeGreen=4',6);