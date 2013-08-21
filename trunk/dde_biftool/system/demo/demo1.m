%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DDE-BIFTOOL demo 1 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add the path to the DDE-BIFTOOL directory, e.g.
%  addpath ../../ddebiftool
disp('Add the correct path ...');
disp('e.g. addpath ../../ddebiftool');
% Note: check whether the path is set in the sys_init() function.

clear;

% init system:

[name,n]=sys_init

%% name = neuron
%% n = 2

% construct a first, approximate steady state point:

stst.kind='stst';
stst.parameter=[1/2 -1 1 2.34 0.2 0.2 1.5];
stst.x=[0 0]'

%% stst = kind: 'stst'
%%   parameter: [0.5000 -1 1 2.3400 0.2000 0.2000 1.5000]
%%           x: [2x1 double]

% get default method parameters for stst calculations:

% we use the new steplength heuristic in the computation of the
% stability of steady states
flag_newhheur=1;

% flag_newhheur=1 is the default choice if this argument is omitted
% from the function df_mthod.
method=df_mthod('stst',flag_newhheur)

method.stability.minimal_real_part=-1;

%% method = continuation: [1x1 struct]
%%                 point: [1x1 struct]
%%             stability: [1x1 struct]

% correct the point:

[stst,success]=p_correc(stst,[],[],method.point)

%% stst = kind: 'stst'
%%   parameter: [0.5000 -1 1 2.3400 0.2000 0.2000 1.5000]
%%           x: [2x1 double]
%% success = 1

stst.x

%% ans = 0
%%       0

% compute its stability:

stst.stability=p_stabil(stst,method.stability);

% plot its stability:

figure(1); clf;
p_splot(stst);

% ask for roots with more negative real part:

method.stability.minimal_real_part=-2;

% recompute stability:

stst.stability=p_stabil(stst,method.stability);

% plot stability:

figure(2); clf;
p_splot(stst);

% get an empty branch with 4 as a free parameter:

branch1=df_brnch(4,'stst')

%% branch1 = method: [1x1 struct]
%%        parameter: [1x1 struct]
%%            point: []

branch1.parameter

%% ans = free: 4
%%  min_bound: [3x2 double]
%%  max_bound: []
%%   max_step: []

branch1.parameter.min_bound 

%% ans = 5  0
%%       6  0
%%       7  0

branch1.parameter.min_bound(4,:)=[4 0];
branch1.parameter.max_bound(1,:)=[4 5];
branch1.parameter.max_step(1,:)=[4 0.2];

% use stst as a first branch point:

branch1.point=stst;

% perturb and correct the point:

stst.parameter(4)=stst.parameter(4)+0.1;
[stst,success]=p_correc(stst,[],[],method.point);

% use as a second branch point:

branch1.point(2)=stst;

% set some continuation parameters:

branch1.method.continuation.plot=0;

% continue in one direction:

[branch1,s,f,r]=br_contn(branch1,100)

%% BR_CONTN warning: boundary hit.
%% branch1 = method: [1x1 struct]
%%        parameter: [1x1 struct]
%%            point: [1x15 struct]
%% s = 15
%% f = 0
%% r = 0

% turn the branch around:

branch1=br_rvers(branch1);

% continue in the other direction:

[branch1,s,f,r]=br_contn(branch1,100);

%% BR_CONTN warning: boundary hit.

% compute stability along branch:

branch1.method.stability.minimal_real_part=-2;
branch1=br_stabl(branch1,0,0);

% obtain suitable scalar measures to plot stability along branch:

[xm,ym]=df_measr(1,branch1);

% plot stability along branch:

figure(3); clf;
br_plot(branch1,xm,ym,'b');
ym

%% ym = field: 'stability'
%%   subfield: 'l1'
%%        row: 'all'
%%        col: 1
%%       func: 'real'

ym.subfield='l0';
br_plot(branch1,xm,ym,'c');
plot([0 5],[0 0],'-.');
axis([0 5 -2 1.5]);

% plot stability versus point number:

figure(4); clf;
br_plot(branch1,[],ym,'b');
br_plot(branch1,[],ym,'b.');
plot([0 30],[0 0],'-.');

% select a point near hopf and turn into hopf guess:

hopf=p_tohopf(branch1.point(24));

% get hopf calculation method parameters:

method=df_mthod('hopf',flag_newhheur);
method.stability.minimal_real_part=-1;

% correct hopf:

[hopf,success]=p_correc(hopf,4,[],method.point)

%% hopf = kind: 'hopf'
%%   parameter: [0.5000 -1 1 0.8071 0.2000 0.2000 1.5000]
%%           x: [2x1 double]
%%           v: [2x1 double]
%%       omega: 0.7820
%% success = 1

% store hopf point in other variable for later use:

first_hopf=hopf;

% compute stability of hopf point:

hopf.stability=p_stabil(hopf,method.stability);

% plot stability of hopf point:

figure(5); clf;
p_splot(hopf);

% use hopf point as first point of hopf branch:

branch2=df_brnch([4 7],'hopf');

branch2.parameter.min_bound(4,:)=[4 0];
branch2.parameter.max_bound(1:2,:)=[[4 4]' [7 10]']';
branch2.parameter.max_step(1:2,:)=[[4 0.2]' [7 0.5]']';

branch2.point=hopf;

% perturb hopf point:

hopf.parameter(7)=hopf.parameter(7)+0.1;

% correct hopf point and recompute stability:

[hopf,success]=p_correc(hopf,4,[],method.point);

% use as second point of hopf branch:

branch2.point(2)=hopf;

% continue with plotting hopf branch:

figure(6); clf;
[branch2,s,f,r]=br_contn(branch2,40);

%% BR_CONTN warning: boundary hit.

branch2=br_rvers(branch2);
[branch2,s,f,r]=br_contn(branch2,20);

% compute stability along branch:

branch2=br_stabl(branch2,0,0);

% plot stability versus point number:

figure(7); clf;
[xm,ym]=df_measr(1,branch2);
ym.subfield='l0';
br_plot(branch2,[],ym,'c');
ym.subfield='l1';
br_plot(branch2,[],ym,'b');

% plot omega to identify 'false' turning point
% as Bogdanov-Takens point:

figure(8); clf;
[xm,ym]=df_measr(0,branch2);
ym

%% ym = field: 'parameter'
%%   subfield: ''
%%        row: 1
%%        col: 7
%%       func: ''

ym.field='omega';
ym.col=1;
xm

%% xm = field: 'parameter'
%%   subfield: ''
%%        row: 1
%%        col: 4
%%       func: ''

xm.col=7;
br_plot(branch2,xm,ym,'c');
grid;

% switch to second hopf near double hopf point:

hopf=p_tohopf(branch2.point(4));

% try to correct second hopf:

[hopf,success]=p_correc(hopf,4,[],method.point)

%% hopf = kind: 'hopf'
%%   parameter: [0.5000 -1 1 -0.0103 0.2000 0.2000 8.5531]
%%           x: [2x1 double]
%%           v: [2x1 double]
%%       omega: 0.9768
%% success = 0

% look at residual behaviour:

method.point.print_residual_info=1;
format short e;
hopf=p_tohopf(branch2.point(4));
[hopf,success]=p_correc(hopf,4,[],method.point);

%% norm_residual = 1.0000e+00  9.3116e-03
%% norm_residual = 2.0000e+00  5.4574e-01
%% norm_residual = 3.0000e+00  6.2629e-02
%% norm_residual = 4.0000e+00  1.8903e-03
%% norm_residual = 5.0000e+00  3.2357e-05

% retry to correct second hopf:

hopf=p_tohopf(branch2.point(4));
[hopf,success]=p_correc(hopf,7,[],method.point)

%% norm_residual = 1.0000e+00   9.3116e-03
%% norm_residual = 2.0000e+00   6.8069e-04
%% norm_residual = 3.0000e+00   2.3169e-07
%% norm_residual = 4.0000e+00   4.3057e-13
%% hopf = kind: 'hopf'
%%   parameter: [5.0000e-01 -1 1 2.0657e-01 2.0000e-01 2.0000e-01 8.6340e+00]
%%           x: [2x1 double]
%%           v: [2x1 double]
%%       omega: 9.1581e-01
%% success = 1

% use as first branch point:

branch3=df_brnch([4 7],'hopf');
branch3.parameter=branch2.parameter;

branch3.point=hopf;

% perturb and correct:

hopf.parameter(4)=hopf.parameter(4)-0.05;
method.point.print_residual_info=0; format short;
[hopf,success]=p_correc(hopf,7,[],method.point);

% use as second branch point:

branch3.point(2)=hopf;

% continue branch of hopf points on two sides:

branch3.method.continuation.plot_progress=0;
figure(6);
[branch3,s,f,r]=br_contn(branch3,100);
%% BR_CONTN warning: boundary hit.
branch3=br_rvers(branch3);
[branch3,s,f,r]=br_contn(branch3,100);
%% BR_CONTN warning: boundary hit.

% turn first hopf point into periodic solution 
% guess with small amplitude 1e-2:

intervals=18;
degree=3;

[psol,stepcond]=p_topsol(first_hopf,1e-2,degree,intervals);

% correct periodic solution guess:

method=df_mthod('psol');
[psol,success]=p_correc(psol,4,stepcond,method.point)

%% psol = kind: 'psol'
%%   parameter: [0.5000 -1 1 0.8072 0.2000 0.2000 1.5000]
%%        mesh: [1x55 double]
%%      degree: 3
%%     profile: [2x55 double]
%%      period: 8.0354
%% success = 1

% empty branch:

branch4=df_brnch(4,'psol');
branch4.parameter.min_bound(4,:)=[4 0];
branch4.parameter.max_bound(1,:)=[4 5];
branch4.parameter.max_step(1,:)=[4 0.1];

% make degenerate periodic solution with amplitude zero
% at hopf point:

deg_psol=p_topsol(first_hopf,0,degree,intervals);

% use these as first two points on branch:
deg_psol.mesh=[];
branch4.point=deg_psol;
psol.mesh=[];
branch4.point(2)=psol;

% compute periodic solutions branch:

figure(9); clf;
[branch4,s,f,r]=br_contn(branch4,50);
%axis([2.3 2.4 0.95 1.15]);

% look at last solution profiles

ll=length(branch4.point);
figure(10); clf;
subplot(3,1,1);
p_pplot(branch4.point(ll-10));
subplot(3,1,2);
p_pplot(branch4.point(ll-5));
subplot(3,1,3);
p_pplot(branch4.point(ll-1));

% look at the period along the branch:

figure(11); clf;
[xm,ym]=df_measr(0,branch4);
ym

%% ym = field: 'profile'
%%   subfield: ''
%%        row: 1
%%        col: 'ampl'
%%       func: ''

ym.field='period';
ym.col=1;
br_plot(branch4,xm,ym,'b');
axis([2.2 2.36 20 170]);

% compute multipliers for point before first turning point:

psol=branch4.point(ll-11);
psol.stability=p_stabil(psol,method.stability);

% plot stability of this point and a point a little further
% (beyond the turning point):

figure(12); clf;
subplot(2,1,1);
p_splot(psol);
axis image;
psol=branch4.point(ll-8);
psol.stability=p_stabil(psol,method.stability);
subplot(2,1,2);
p_splot(psol);

% recompute on a better mesh :

psol=branch4.point(ll-12);
intervals=40;
degree=4;
psol=p_remesh(psol,degree,intervals);
method.point.adapt_mesh_after_correct=1;
method.point.newton_max_iterations=7;
method.point.newton_nmon_iterations=2;
[psol,success]=p_correc(psol,[],[],method.point)

%% psol = kind: 'psol'
%%   parameter: [0.5000 -1 1 2.3358 0.2000 0.2000 1.5000]
%%        mesh: [1x161 double]
%%      degree: 4
%%     profile: [2x161 double]
%%      period: 38.4916
%% success = 1

% use this point to start a branch:

branch5=df_brnch(4,'psol');
branch5.parameter=branch4.parameter;

branch5.point=psol;
psol.parameter(4)=psol.parameter(4)+0.01;
[psol,success]=p_correc(psol,[],[],method.point,1);
branch5.point(2)=psol;
branch5.method=method;
[xm,ym]=df_measr(0,branch5);
ym.field='period';
ym.col=1;
figure(11); axis auto; hold on;
branch5.method.continuation.plot_measure.x=xm;
branch5.method.continuation.plot_measure.y=ym;
[branch5,s,f,r]=br_contn(branch5,25);

% stability beyond the turning point:

psol=branch5.point(6);
psol.stability=p_stabil(psol,method.stability);
psol.stability.mu

%% ans = 241.2300
%%         1.0000

% plotting the profile reveals the double homoclinic structure:

figure(13); clf;
subplot(2,1,1);
ll=length(branch5.point);
psol=branch5.point(ll-5);
plot(psol.mesh,psol.profile);
subplot(2,1,2);
psol1=p_remesh(psol,degree,0:0.001:1);
psol2=p_remesh(psol,degree,(0:0.001:1)+0.02);
plot(psol1.profile',psol2.profile');
psol.period

%% ans = 399.7466


% we take the first half of the profile and rescale to [0,1]
figure(14);clf;subplot(2,1,1);
hcli1=psol;
hcli1.mesh=hcli1.mesh(1:65);
hcli1.profile=hcli1.profile(:,1:65);
hcli1.period=hcli1.period*hcli1.mesh(end);
hcli1.mesh=hcli1.mesh/hcli1.mesh(end);

% convert it to a point of homoclinic structure
hcli1=p_tohcli(hcli1)

%% hcli1 = 
%%         kind: 'hcli'
%%    parameter: [0.5000 -1 1 2.3460 0.2000 0.2000 1.5000]
%%         mesh: [1x61 double]
%%       degree: 4
%%      profile: [2x61 double]
%%       period: 113.4318
%%           x1: [2x1 double]
%%           x2: [2x1 double]
%%     lambda_v: 0.3142
%%     lambda_w: 0.3142
%%            v: [2x1 double]
%%            w: [2x1 double]
%%        alpha: 1
%%      epsilon: 2.9010e-04

% and correct it
mh=df_mthod('hcli');
[hcli1,success]=p_correc(hcli1,4,[],mh.point)

%% hcli1 = 
%%         kind: 'hcli'
%%    parameter: [0.5000 -1 1 2.3460 0.2000 0.2000 1.5000]
%%         mesh: [1x61 double]
%%       degree: 4
%%      profile: [2x61 double]
%%       period: 114.8378
%%           x1: [2x1 double]
%%           x2: [2x1 double]
%%     lambda_v: 0.3141
%%     lambda_w: 0.3141
%%            v: [2x1 double]
%%            w: [2x1 double]
%%        alpha: 1
%%      epsilon: 2.9010e-04
%% success =
%%     1

p_pplot(hcli1);

% we can do the same for the second half of the profile
figure(14);subplot(2,1,2);
hcli2=psol;
hcli2.mesh=hcli2.mesh(81:end-16);
hcli2.profile=hcli2.profile(:,81:end-16);
hcli2.mesh=hcli2.mesh-hcli2.mesh(1);
hcli2.period=hcli2.period*hcli2.mesh(end);
hcli2.mesh=hcli2.mesh/hcli2.mesh(end);

hcli2=p_tohcli(hcli2);
[hcli2,success]=p_correc(hcli2,4,[],mh.point);
p_pplot(hcli2);

% we recompute hcli1 with more mesh points
hcli1=p_remesh(hcli1,4,70);
[hcli1,success]=p_correc(hcli1,4,[],mh.point)

%% hcli1 = 
%%         kind: 'hcli'
%%    parameter: [0.5000 -1 1 2.3460 0.2000 0.2000 1.5000]
%%         mesh: [1x281 double]
%%       degree: 4
%%      profile: [2x281 double]
%%       period: 115.5581
%%           x1: [2x1 double]
%%           x2: [2x1 double]
%%     lambda_v: 0.3142
%%     lambda_w: 0.3142
%%            v: [2x1 double]
%%            w: [2x1 double]
%%        alpha: 1
%%      epsilon: 2.9010e-04
%% success =
%%     1

% we continue the the first homoclinic orbit with respect to tau_s
figure(15);
branch6=df_brnch([4 7],'hcli');
branch6.point=hcli1;
hcli1.parameter(7)=1.49;
[hcli1,success]=p_correc(hcli1,4,[],mh.point);
branch6.point(2)=hcli1;
[branch6,s,r,f]=br_contn(branch6,19);

% save data:

save demo1 branch1 branch2 branch3 branch4 branch5 branch6;

% exit;
