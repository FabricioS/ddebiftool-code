clear;                           % clear variables
close all;                       % close figures
addpath('../../../ddebiftool',...
    '../../../ddebiftool_extra_psol',...
    '../../../ddebiftool_extra_nmfm',...
    '../../../ddebiftool_utilities');

funcs=set_funcs(...
    'sys_rhs', @(xx,par) sys_rhs(xx,par),...
    'sys_tau', @() 2,...
    'sys_deri', @(xx,par,nx,np,v) sys_deri(xx,par,nx,np,v),...
    'sys_mfderi',@(xx,par,varargin) sys_mfderi(xx,par,varargin{:}));

disp('Holling-Tanner');
% par=(beta,tau,a,m,h,delta)

%% staility of a BT-point
beta=0.5;alpha=0.5;
m=(1/30)*(1-beta/(alpha*beta+1));
h=(1/4)*(beta/(alpha*beta+1)-1)^2+m*beta/(alpha*beta+1);
tau=1/4*(alpha*beta+1)^2/beta;
delta=1/(alpha*beta+1)^2;

xster=-(1/2)*(beta/(alpha*beta+1)+2*m-1);
yster=beta*xster;

p.parameter=[beta,tau,alpha,m,h,delta];
p.x=[xster;yster];

p.kind='stst';

method=df_mthod(funcs,'stst');
method.stability.minimal_real_part=-1;

% disp('before');
% disp(p.x)

% [p,success]=p_correc(funcs,p,[],[],method.point);

p.stability=p_stabil(funcs,p,method.stability);
p.stability.l1

% [p,success]=p_correc(funcs,p,[],[],method.point);


%% Continuation of steady state branch
beta=0.5;alpha=0.5;
m=(1/30)*(1-beta/(alpha*beta+1));
h=(1/4)*(beta/(alpha*beta+1)-1)^2+m*beta/(alpha*beta+1);
tau=1/4*(alpha*beta+1)^2/beta;
delta=1/(alpha*beta+1)^2;

beta=0.486844;
delta=0.5;
stst.parameter=[beta,tau,alpha,m,h,delta]

% xster=-(1/2)*(beta/(alpha*beta+1)+2*m-1)
xster=0.3331;
yster=beta*xster;

contpar=1;
stst_branch = SetupStst(funcs,'x',[xster; yster],'parameter',stst.parameter,...
    'contpar',contpar,'max_step',[contpar 0.002],'min_bound',...
    [contpar 0.4],'max_bound',[contpar 0.6],...
    'newheuristics_tests',0);

stst_branch.method.continuation.plot = 1;

[stst_branch,s,f,r] = br_contn(funcs,stst_branch,100);
stst_branch = br_rvers(stst_branch);
[stst_branch,s,f,r] = br_contn(funcs,stst_branch,100);

stst_branch=br_bifdet(funcs,stst_branch);

%% stability of Hopf point
point_nr=10;
p.x = stst_branch.point(point_nr).x;
p.parameter = stst_branch.point(point_nr).parameter;
p.kind='stst';
method=df_mthod(funcs,'stst');
% method.stability.minimal_real_part=-1;

p.stability=p_stabil(funcs,p,method.stability);
p.stability.l1


%% plot stability along the stst_branch
stst_branch = br_stabl(funcs,stst_branch,0,0);
[xm,ym] = df_measr(1,stst_branch);
ym.subfield='l0';
ym.func='real';
br_plot(stst_branch,[],ym,'-');

%% plot stability along the stst_branch
stst_branch = br_stabl(funcs,stst_branch,0,0);
[xm,ym] = df_measr(0,stst_branch);
br_plot(stst_branch,xm,[],'-');

%% Continuation of the first Hopf branch
FPI = br_getflags(stst_branch);
start_ind = FPI(1);
% start_ind = ind_hopf;
% start_ind=54;

% We do a standard continuation, using the starting index obtained from the
% flagged point indices.
fprintf('----- Hopf branch 1 -----\n');
[hopf_branch, suc] = SetupHopf(funcs, stst_branch, start_ind, 'contpar', [1 6], 'dir', 1, 'step', 0.002);

betamin=0.4; betamax=0.6; deltamin=0.4; deltamax=0.7;

hopf_branch.parameter.min_bound=[1 betamin; 6 deltamin];
hopf_branch.parameter.max_bound=[1 betamax; 6 deltamax];
hopf_branch.parameter.max_step=[1 0.005; 6 0.005];

hopf_branch.method.continuation.plot = 0;

[hopf_branch,s,f,r]=br_contn(funcs,hopf_branch,250);

[hopf_branch,success] = br_bifdet(funcs,hopf_branch);


%% Continuation plot of Hopf branch
close all;
hopf_branch  = br_stabl(funcs,hopf_branch,0,0);
[xm,ym] = df_measr(0,hopf_branch);

% default
br_plot(hopf_branch,xm,ym,'-*')
xlabel('$\beta$','interpreter','latex');
ylabel('$\delta$','interpreter','latex');

% point number vs delta
figure;
br_plot(hopf_branch,[],ym,'-*')
xlabel('point number');
ylabel('$\delta$','interpreter','latex');


%% get periodic orbit emanating from the hopf curve
hopf=hopf_branch.point(10);
[psol,stp]=p_topsol(funcs,hopf,1e-2,3,27);

ind_delta=6;
mpsol=df_mthod(funcs,'psol');
[psol,s]=p_correc(funcs,psol,ind_delta,stp,mpsol.point);
psol_branch=df_brnch(funcs,ind_delta,'psol');
psol_branch.point=psol;
[psol,stp]=p_topsol(funcs,hopf,1.1e-2,3,27);
[psol,s]=p_correc(funcs,psol,ind_delta,stp,mpsol.point);
psol_branch.point(2)=psol;
figure(2);clf;
[xm,ym]=df_measr(0,psol_branch);
ym.field='period';
ym.col=1;
ym.row=1;
psol_branch.method.continuation.plot_measure.x=xm;
psol_branch.method.continuation.plot_measure.y=ym;
[psol_branch,s,r,f]=br_contn(funcs,psol_branch,50);
xlabel('$\delta$','interpreter','latex');
ylabel('period');

%%
figure(3);clf;
psol=psol_branch.point(end)
p_pplot(psol);
xlabel('time/period');ylabel('x1,x2');

%% convert last point of the psol_branch to homoclinic orbit
hcli=p_tohcli(funcs,psol)
figure(4);clf;
p_pplot(hcli);
xlabel('time/period');ylabel('x1,x2');

%% correct homoclinic orbit
mhcli=df_mthod(funcs,'hcli');
[hcli,s]=p_correc(funcs,hcli,ind_delta,[],mhcli.point); % correct
hcli=p_remesh(hcli,3,50); % remesh it and
[hcli,s]=p_correc(funcs,hcli,ind_delta,[],mhcli.point); % correct it again
figure(5);clf;
p_pplot(hcli) % plot it after remeshing & correction
xlabel('time/period');ylabel('x1,x2');

%% perturb hcli, correct and continue
ind_beta=1;
hcli_br=df_brnch(funcs,[ind_beta, ind_delta],'hcli');
hcli_br.point=hcli;
hcli.parameter(ind_beta)=hcli.parameter(ind_beta)-6e-6;
[hcli,s]=p_correc(funcs,hcli,ind_delta,[],mhcli.point);
hcli_br.point(2)=hcli;

figure(1);
[hcli_br,s,r,f]=br_contn(funcs,hcli_br,28);
hcli_br=br_rvers(hcli_br);
[hcli_br,s,r,f]=br_contn(funcs,hcli_br,40);

xlabel('$\beta$','interpreter','latex');
ylabel('$\delta$','interpreter','latex');
%% BT-point
figure(8);
stst=p_tostst(funcs,hcli_br.point(end));
stst=stst(1);
mstst=df_mthod(funcs,'stst');
stst.stability=p_stabil(funcs,stst,mstst.stability);
p_splot(stst);
xlabel('\Re\lambda');ylabel('\Im\lambda');

%% continue fold
free_par=1;
fold1=p_tofold(funcs,stst_branch.point(49));
[fold1,success]=p_correc(funcs,fold1,free_par,[],method.point);

fold2 = p_axpy(0.9999,fold1,[]); % perturb fold point to get a secondary point for the branch

fold2=fold1;
fold2.parameter(1)=0.5;
[fold2,success]=p_correc(funcs,fold2,free_par,[],method.point); % correct fold point
fold2.parameter(1)=0.5;
% fold2.stability = p_stabil(funcs,fold2,method.stability);

fold_branch=df_brnch(funcs,free_par,'fold'); % use fold point as first point of fold branch:
fold_branch.point=fold1;
fold_branch.point(2)=fold2;

fold_branch.parameter.max_step(1:2,:)=[[1 0.01]' [6 0.01]']';

fold_branch.method.continuation.plot = 1;

[fold_branch,s,f,r]=br_contn(funcs,fold_branch,10);
fold_branch=br_rvers(fold_branch);
[fold_branch,s,f,r]=br_contn(funcs,fold_branch,100);

