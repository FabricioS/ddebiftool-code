%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DDE-BIFTOOL hom_demo %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
format short e;
sys_init;

% loading and plotting branches of folds and hopf points
load hom_demo;
figure(1);
[xm,ym]=df_measr(0,fold_branch);
br_plot(fold_branch,xm,ym,':');
axis([1.28 1.62 -1.36 -1.24]);
hold on;
br_plot(hopf1_branch,xm,ym,'-.');
br_plot(hopf2_branch,xm,ym,'-.');

% select a Hopf point somewhere on one branch, and start the branch of
% periodic solutions that emanates from it.
hopf=hopf1_branch.point(27);
[psol,stp]=p_topsol(hopf,1e-2,3,27)

%% psol = kind: 'psol'
%%   parameter: [2.6000e+00 1.3428e+00 1 -1.3398e+00 -5.0000e-01 1]
%%        mesh: [1x82 double]
%%      degree: 3
%%     profile: [2x82 double]
%%      period: 5.5271e+01
%% stp = kind: 'psol'
%%  parameter: [0 0 0 0 0 0]
%%       mesh: [1x82 double]
%%     degree: 3
%%    profile: [2x82 double]
%%     period: 0

mpsol=df_mthod('psol');
[psol,s]=p_correc(psol,4,stp,mpsol.point);
psol_branch=df_brnch(4,'psol');
psol_branch.point=psol;
[psol,stp]=p_topsol(hopf,2e-2,3,27);
[psol,s]=p_correc(psol,4,stp,mpsol.point);
psol_branch.point(2)=psol;
figure(2);clf;
[xm,ym]=df_measr(0,psol_branch);
ym.field='period';
ym.col=1;
ym.row=1;
psol_branch.method.continuation.plot_measure.x=xm;
psol_branch.method.continuation.plot_measure.y=ym;
[psol_branch,s,r,f]=br_contn(psol_branch,20);

% the last point is close to a homoclinic, as this plot shows.
figure(3);clf;
psol=psol_branch.point(end)
%% psol = kind: 'psol'
%%   parameter: [2.6000e+00 1.3428e+00 1 -1.3392e+00 -5.0000e-01 1]
%%        mesh: [1x82 double]
%%      degree: 3
%%     profile: [2x82 double]
%%      period: 2.1469e+02

p_pplot(psol);

% so we convert it to a point of homoclinic type.
hcli=p_tohcli(psol)

%% hcli = kind: 'hcli'
%%   parameter: [2.6000e+00 1.3428e+00 1 -1.3392e+00 -5.0000e-01 1]
%%        mesh: [1x79 double]
%%      degree: 3
%%     profile: [2x79 double]
%%      period: 1.8216e+02
%%          x1: [2x1 double]
%%          x2: [2x1 double]
%%    lambda_v: 1.6906e-01
%%    lambda_w: 1.6906e-01
%%           v: [2x1 double]
%%           w: [2x1 double]
%%       alpha: -1
%%     epsilon: 5.2583e-06

% we plot it before correct
figure(4);clf;
p_pplot(hcli);

% correct it
mhcli=df_mthod('hcli');
[hcli,s]=p_correc(hcli,4,[],mhcli.point);
% remesh it and correct it
hcli=p_remesh(hcli,3,50);
[hcli,s]=p_correc(hcli,4,[],mhcli.point)

%% hcli = kind: 'hcli'
%%   parameter: [2.6000e+00 1.3428e+00 1 -1.3392e+00 -5.0000e-01 1]
%%        mesh: [1x151 double]
%%      degree: 3
%%     profile: [2x151 double]
%%      period: 1.8806e+02
%%          x1: [2x1 double]
%%          x2: [2x1 double]
%%    lambda_v: 1.6905e-01
%%    lambda_w: 1.6905e-01
%%           v: [2x1 double]
%%           w: [2x1 double]
%%       alpha: -1
%%     epsilon: 5.2583e-06
%% s =
%%      1
 
% and plot it after correct
figure(5);clf;
p_pplot(hcli);

% we now continue a branch of homoclinic solutions in two-parameter
% space, and show this on the first figure.
hcli_br=df_brnch([2 4],'hcli');
hcli_br.point=hcli;
hcli.parameter(2)=hcli.parameter(2)+6e-3;
[hcli,s]=p_correc(hcli,4,[],mhcli.point);
hcli_br.point(2)=hcli;
figure(1);
[hcli_br,s,r,f]=br_contn(hcli_br,85)
%% hcli_br = 
%%        method: [1x1 struct]
%%     parameter: [1x1 struct]
%%         point: [1x71 struct]
%% s = 70
%% r = 16
%% f = 0
hcli_br=br_rvers(hcli_br);
[hcli_br,s,r,f]=br_contn(hcli_br,10)
%% hcli_br = 
%%        method: [1x1 struct]
%%     parameter: [1x1 struct]
%%         point: [1x81 struct]
%% s = 11
%% r = 0
%% f = 0

% we now do the same for the second branch
hopf=hopf2_branch.point(27);
[psol,stp]=p_topsol(hopf,1e-2,3,27);
mpsol=df_mthod('psol');
[psol,s]=p_correc(psol,4,stp,mpsol.point);
psol_branch=df_brnch(4,'psol');
psol_branch.point=psol;
[psol,stp]=p_topsol(hopf,2e-2,3,27);
[psol,s]=p_correc(psol,4,stp,mpsol.point);
psol_branch.point(2)=psol;
psol_branch.method.continuation.plot=0;
psol_branch.method.continuation.plot_progress=0;
[psol_branch,s,r,f]=br_contn(psol_branch,20);

% again the last point is close to a homoclinic
% so we convert it to a point of homoclinic type.

psol=psol_branch.point(end);
hcli=p_tohcli(psol);

% correct it
mhcli=df_mthod('hcli');
[hcli,s]=p_correc(hcli,4,[],mhcli.point);

% remesh it and correct it
hcli=p_remesh(hcli,3,50);
[hcli,s]=p_correc(hcli,4,[],mhcli.point);

% we now continue this second  branch of homoclinic solutions in two-parameter
% space, and show this on the first figure.
hcli2_br=df_brnch([2 4],'hcli');
hcli2_br.point=hcli;
hcli.parameter(2)=hcli.parameter(2)+6e-3;
[hcli,s]=p_correc(hcli,4,[],mhcli.point);
hcli2_br.point(2)=hcli;
figure(1);
[hcli2_br,s,r,f]=br_contn(hcli2_br,85)
%% hcli2_br = 
%%        method: [1x1 struct]
%%     parameter: [1x1 struct]
%%         point: [1x70 struct]
%% s = 69
%% r = 17
%% f = 0
hcli2_br=br_rvers(hcli2_br);
[hcli2_br,s,r,f]=br_contn(hcli2_br,10)
%% hcli2_br = 
%%        method: [1x1 struct]
%%     parameter: [1x1 struct]
%%         point: [1x80 struct]
%% s = 11
%% r = 0
%% f = 0

% we show that for e1=-1.3, we have a double homoclinic orbit
figure(6);
[xm,ym]=df_measr(0,hcli_br);
hold on;
br_plot(hcli2_br,[],ym);
br_plot(hcli_br,[],ym,'--');
plot([0 110],[-1.3 -1.3],'r-.');
axis([22 40 -1.304 -1.294]);

figure(7);
plot(hcli2_br.point(30).profile(1,:),hcli2_br.point(30).profile(2,:));
hold on;
plot(hcli_br.point(26).profile(1,:),hcli_br.point(26).profile(2,:));

% each branch of homoclinic orbits emanates from a Takens-Bogdanov point.
% Such a point has a double zero eigenvalue.  
figure(8);
stst=p_tostst(hcli_br.point(end));
stst=stst(1);
mstst=df_mthod('stst');
stst.stability=p_stabil(stst,mstst.stability);
p_splot(stst);

%try to get closer by initializing a new branch and taking very small steps
hcli=hcli_br.point(end);
hcli=p_remesh(hcli,3,70);
[hcli,s]=p_correc(hcli,4,[],mhcli.point);
to_tb_branch=df_brnch([2 4],'hcli');
to_tb_branch.point=hcli;
hcli.parameter(2)=hcli.parameter(2)-1e-4;
hcli=p_remesh(hcli,3,70);
[hcli,s]=p_correc(hcli,4,[],mhcli.point);
to_tb_branch.point(2)=hcli;

to_tb_branch.method.continuation.plot=0;
to_tb_branch.method.continuation.plot_progress=0;
[to_tb_branch,s,r,f]=br_contn(to_tb_branch,40);

figure(9);
stst=p_tostst(to_tb_branch.point(end));
stst=stst(1);
mstst=df_mthod('stst');
stst.stability=p_stabil(stst,mstst.stability);
p_splot(stst);

save hom_demo;

