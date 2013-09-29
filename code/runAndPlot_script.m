clc; clear all; close all;
randn('state',123); rand('state',0);

params.n = 300;
params.p = 40;
params.noise = 0.05;
params.s = 3;
params.gamma = 0.5;

params.fig_name = 'QuadLin3';
params.use_quad = 1;
params.use_lin = 1;
params.use_cube = 0;
params.use_logexp = 0;

out = runAndPlot(params);


%==========%

clc; clear all; close all;
randn('state',123); rand('state',0);

params.n = 120;
params.p = 40;
params.noise = 0.05;
params.s = 3;
params.gamma = 0.05;

params.fig_name = 'QuadLin2';
params.use_quad = 1;
params.use_lin = 1;
params.use_cube = 0;
params.use_logexp = 0;

out = runAndPlot(params);

%==========%

clc; clear all; close all;
randn('state',123); rand('state',0);

params.n = 120;
params.p = 40;
params.noise = 0.05;
params.s = 3;
params.gamma = 0.5;

params.fig_name = 'QuadLin';
params.use_quad = 1;
params.use_lin = 1;
params.use_cube = 0;
params.use_logexp = 0;

out = runAndPlot(params);

%==========%

clc; clear all; close all;
randn('state',123); rand('state',0);

params.n = 100;
params.p = 40;
params.noise = 0.02;
params.s = 3;
params.gamma = 0.5;

params.fig_name = 'Quadratic';
params.use_quad = 1;
params.use_lin = 0;
params.use_cube = 0;
params.use_logexp = 0;

out = runAndPlot(params);

%==========%

clc; clear all; close all;
randn('state',123); rand('state',0);

params.n = 100;
params.p = 40;
params.noise = 0.02;
params.s = 3;
params.gamma = 0.5;

params.fig_name = 'Logexp';
params.use_quad = 0;
params.use_lin = 0;
params.use_cube = 0;
params.use_logexp = 1;

out = runAndPlot(params);

%==========%

clc; clear all; close all;
randn('state',123); rand('state',0);

params.n = 100;
params.p = 40;
params.noise = 0.02;
params.s = 3;
params.gamma = 0.5;

params.fig_name = 'Lin';
params.use_quad = 0;
params.use_lin = 1;
params.use_cube = 0;
params.use_logexp = 0;

out = runAndPlot(params);

%==========%

clc; clear all; close all;
randn('state',123); rand('state',0);

params.n = 100;
params.p = 40;
params.noise = 0.02;
params.s = 3;
params.gamma = 0.5;

params.fig_name = 'Cube';
params.use_quad = 0;
params.use_lin = 0;
params.use_cube = 1;
params.use_logexp = 0;

out = runAndPlot(params);