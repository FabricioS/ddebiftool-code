%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DDE-BIFTOOL sd-demo %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

% init system:

[name,n]=sys_init

% construct a first, approximate steady state point:

stst.kind='stst';
stst.parameter=[4.5 0.04 -1.4 6 -0.45 -0.01 3 0.3 0.1 1 0.2];
stst.x=[1.4 1.5 -25 0.6 1.4]';

% get default method parameters for stst calculations:

method=df_mthod('stst');

% correct the point:
[stst,success]=p_correc(stst,[],[],method.point)

%% stst = kind: 'stst'
%%   parameter: [1x11 double]
%%           x: [5x1 double]
%% success = 1

stst.x

%%ans = 1.4134
%%      1.5193
%%    -25.1077
%%      0.5886
%%      1.3801

% get an empty branch with 5 as a free parameter:

branch1=df_brnch(5,'stst');

branch1.parameter.min_bound(1,:)=[5 -1];
branch1.parameter.max_bound(1,:)=[5 1];
branch1.parameter.max_step(1,:)=[5 0.1];    

% use stst as a first branch point:

branch1.point=stst;

% perturb and correct the point:

stst.parameter(5)=stst.parameter(5)-0.01;
[stst,success]=p_correc(stst,[],[],method.point);

% use as a second branch point:

branch1.point(2)=stst;

% continue with plotting the branch:

figure(1); clf;
[branch1,s,f,r]=br_contn(branch1,20)

%% BR_CONTN warning: delay number_3 becomes negative.
%% branch1 = method: [1x1 struct]
%%        parameter: [1x1 struct]
%%            point: [1x9 struct]
%% s = 8
%% f = 0
%% r = 0

plot(branch1.point(end).parameter(5),branch1.point(end).x(1),'o');

% look at the value of delay 3 at the last computed point

p_tau(branch1.point(end),3)

%% ans = 2.2204e-16

% compute stability along branch:

branch1.method.stability.minimal_real_part=-1;
branch1=br_stabl(branch1,0,0);

% obtain suitable scalar measures to plot stability along branch:

[xm,ym]=df_measr(1,branch1);
ym.subfield='l1';

% plot stability along branch versus point number:

figure(2); clf;
br_plot(branch1,[],ym,'b');
br_plot(branch1,[],ym,'b.');
plot([0 10],[0 0],'-.');

branch1.point(5).stability.l1

%% ans = -0.0023 - 0.5488i
%%       -0.0023 + 0.5488i
%%       -0.0952          
%%       -0.4499  

% select a point near hopf and turn into hopf guess:

hopf=p_tohopf(branch1.point(5));

% get hopf calculation method parameters:

method=df_mthod('hopf');

% correct hopf:

[hopf,success]=p_correc(hopf,5,[],method.point)

%% hopf = kind: 'hopf'
%%   parameter: [1x11 double]
%%           x: [5x1 double]
%%           v: [5x1 double]
%%       omega: 0.5497
%% success = 1

% use hopf point as first point of hopf branch:

branch2=df_brnch([2 9],'hopf');

branch2.parameter.min_bound(1:2,:)=[[2 -1]' [9 -1]']';
branch2.parameter.max_bound(1:2,:)=[[2 10]' [9 10]']';
branch2.parameter.max_step(1:2,:)=[[2 1]' [9 1]']';    

branch2.point=hopf;

% perturb hopf point:

hopf.parameter(9)=hopf.parameter(9)+0.1;

% correct hopf point:

[hopf,success]=p_correc(hopf,2,[],method.point);

% use as second point of hopf branch:

branch2.point(2)=hopf;

% continue with plotting hopf branch:

figure(3); clf;
[branch2,s,f,r]=br_contn(branch2,14);

% turn first hopf point into periodic solution 
% guess with small amplitude 0.1:

hopf=branch2.point(1);

intervals=15;
degree=3;

[psol,stepcond]=p_topsol(hopf,1e-1,degree,intervals); 

% correct periodic solution guess:

method=df_mthod('psol');
[psol,success]=p_correc(psol,10,stepcond,method.point)

%% psol = kind: 'psol'
%%   parameter: [1x11 double]
%%        mesh: [1x46 double]
%%      degree: 3
%%     profile: [5x46 double]
%%      period: 11.4306
%% success = 1

% empty branch:
branch3=df_brnch(10,'psol');
branch3.parameter.min_bound(1,:)=[10 0];
branch3.parameter.max_bound(1,:)=[10 10];
branch3.parameter.max_step(1,:)=[10 0.01];

deg_psol=p_topsol(hopf,0,degree,intervals);

% use these as first two points on branch:

branch3.point=deg_psol;
branch3.point(2)=psol;

% compute periodic solutions branch:

figure(4); clf;
[branch3,s,f,r]=br_contn(branch3,10); 

%% BR_CONTN warning: delay number_3 becomes negative.

% indicate the last computed point on the branch
point=branch3.point(end);
p_ampl=max(point.profile(1,:))-min(point.profile(1,:));
plot(point.parameter(10),p_ampl,'o');

%% BR_CONTN warning: delay number_3 becomes negative.

% compute and plot delay 3 at representation points at the last point

tau_eva=p_tau(branch3.point(end),3);
figure(5); clf;
plot(branch3.point(end).mesh,tau_eva);
hold;
plot(branch3.point(end).mesh,tau_eva,'.');

% compute minimum delay 3
min(tau_eva) 

%% ans = 9.6557e-04

% turn last hopf point into periodic solution 
% guess with small amplitude 0.1:

hopf=branch2.point(end);

intervals=15;
degree=3;

[psol,stepcond]=p_topsol(hopf,0.1,degree,intervals);

% correct periodic solution guess:

method=df_mthod('psol');
[psol,success]=p_correc(psol,1,stepcond,method.point)

%% psol = kind: 'psol'
%%   parameter: [1x11 double]
%%        mesh: [1x46 double]
%%      degree: 3
%%     profile: [5x46 double]
%%      period: 12.6610
%% success = 1

% empty branch:
branch4=df_brnch(1,'psol');
branch4.parameter.min_bound(1,:)=[1 0];
branch4.parameter.max_bound(1,:)=[1 10];
branch4.parameter.max_step(1,:)=[1 0.01];

deg_psol=p_topsol(hopf,0,degree,intervals);

% use these as first two points on branch:
branch4.point=deg_psol;
branch4.point(2)=psol;

% compute periodic solutions branch:

figure(6); clf;
[branch4,s,f,r]=br_contn(branch4,20);

%% BR_CONTN warning: delay number_6 becomes negative.

% indicate the last computed point on the branch

point=branch4.point(end);
p_ampl=max(point.profile(1,:))-min(point.profile(1,:));
plot(point.parameter(1),p_ampl,'o');

% plot delay 6 at representation points at the last point

psol=branch4.point(end);
figure(7); clf;
plot(psol.mesh,psol.profile(5,:));
hold;
plot(psol.mesh,psol.profile(5,:),'.');
min(psol.profile(5,:)) 

%% ans = -5.8556e-31

% compute and plot stability of solution at the last point

psol.stability=p_stabil(psol,method.stability);
psol.stability.mu
%% ans = 1.3253
%%       1.0000
%%       0.0959

figure(8); clf;
p_splot(psol);
axis image;

% save data:

save sd_demo branch1 branch2 branch3 branch4;

% exit;

