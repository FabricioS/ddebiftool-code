function hcli=p_tohcli(point)

% INPUT:
%     point a periodic solution near a homoclinic solution
%           alternatively an initial point in a hcli structure,
%           where a good starting guess for the profile and steady
%           states are available
% OUTPUT:
%     hcli a starting value to compute the exact homoclinic or
%     heteroclinic solution  

% (c) DDE-BIFTOOL v. 2.02, 16/6/2002

    if mod(length(point.mesh),point.degree)~=1,
      err=[length(point.mesh) point.degree];
      error('P_TOHCLI: psol does not contain L intervals of m points!');
    end;
    
    hcli.kind='hcli';
    hcli.parameter=point.parameter;
    hcli.mesh=point.mesh;
    hcli.degree=point.degree;
    
    switch point.kind,
     
     case 'psol',
      
      for i=1:length(point.profile(1,:))-1,
	test(i)=norm(point.profile(:,i)-point.profile(:,i+1));
      end;
      [minval pos]=min(abs(test));
      stst.kind='stst';
      stst.parameter=hcli.parameter;
      stst.x=point.profile(:,pos);
      
      for i=1:size(point.profile,2)
	x_profile(1,i)=norm(point.profile(:,i)-stst.x);
      end;
      
      [peak peak_pos]=max(x_profile(1,:));
      [hole hole_pos]=min(x_profile(1,:));
      left_part=point.profile(:,1:peak_pos);
      right_part=point.profile(:,peak_pos+1:end);
      hole_begin=hole_pos-mod(hole_pos,point.degree)+1;
      hole_end=hole_begin+point.degree;
      
      if hole_pos<peak_pos,
	right_part=[right_part left_part(:,2:hole_begin)];
	left_part=left_part(:,hole_end:end);
	hcli.mesh=[hcli.mesh(hole_end:end) ...
		   (hcli.mesh(2:hole_begin)+1)];
      else
	left_part=[right_part(:,hole_end-peak_pos:end-1) left_part];
	right_part=right_part(:,1:hole_begin-peak_pos);
	hcli.mesh= ...
	    [(hcli.mesh(hole_end:end-1)-1)...
	     hcli.mesh(1:hole_begin)];
      end;
      
      nb_of_points=length(hcli.mesh);
      rest=mod(nb_of_points,point.degree);
      hcli.profile=[left_part right_part];
      
      if rest>1,
	hcli.profile=point.profile(:,1+floor((rest-1)/2):end-ceil((rest-1)/2));
	hcli.mesh=hcli.mesh(1+floor((rest-1)/2):end-ceil((rest-1)/2));
      end;
      if rest==0,
	rest=point.degree;
	hcli.profile=point.profile(:,1+floor((rest-1)/2):end-ceil((rest-1)/2));
	hcli.mesh=hcli.mesh(1+floor((rest-1)/2):end-ceil((rest-1)/2));
      end;
      hcli.mesh=hcli.mesh-hcli.mesh(1);
      hcli.period=point.period*hcli.mesh(end);
      hcli.mesh=hcli.mesh/hcli.mesh(end);
      hcli.x1=stst.x;
      hcli.x2=stst.x;
      stst1=stst;
      stst2=stst;
     case 'hcli',
      hcli=point;
      stst1.kind='stst';
      stst1.parameter=point.parameter;
      stst1.x=point.x1;
      stst2=stst1;
      stst2.x=point.x2;
     otherwise,
      err=1;
      error(['P_TOHCLI: not a valid conversion for other than psol' ...
	     ' or hcli type points']);
    end;
    
    m=df_mthod('stst');
    stst1.stability=p_stabil(stst1,m.stability);
    
    i=1;
    if stst1.stability.l1(i)<0,
      err=stst.stability.l1;
      error('P_TOHCLI: no unstable eigenmodes found');
    end;
    while (i<= length(stst1.stability.l1) & real(stst1.stability.l1(i))>0 )
      lambda(i,1)=stst1.stability.l1(i); 
      i=i+1;
    end;
    
    hcli.lambda_v=lambda;
    
    tp_del=nargin('sys_tau');
    if tp_del==0
      tau=point.parameter(sys_tau);
      n_tau=length(tau);
    else
      error('P_TOHCLI: computing connected orbits is not implemented for equations with state-dependent delays');
    end;
    
    
    n=length(hcli.profile(:,1));
    if n==1,
      v=ones(1,length(lambda));
    else
      for i=1:length(lambda)
	delta=eye(n)*lambda(i);
	xx=stst1.x;
	for t=1:n_tau
	  xx=[xx stst1.x];
	end;
	delta=delta-sys_deri(xx,hcli.parameter,0,[],[]);
	for t=1:n_tau
	  delta=delta-sys_deri(xx,hcli.parameter,t,[],[])*exp(-lambda(i)*tau(t));
	end;
	[eigvec,eigval]=eig(delta);
	[minval pos]=min(abs(diag(eigval)));
	v(:,i)=eigvec(:,pos);
      end;
    end;
    
   stst2.stability=p_stabil(stst2,m.stability);
   
   lambda=[];
    
    i=1;
    while (i<= length(stst2.stability.l1) & real(stst2.stability.l1(i))>0 )
      lambda(i,1)=stst2.stability.l1(i); 
      i=i+1;
    end;
   
    hcli.lambda_w=lambda;
    
    if n==1
      w=ones(1,length(lambda));
    else
      for i=1:length(lambda),
	     delta=eye(n)*lambda(i);
	     xx=stst2.x;
	     for t=1:n_tau
	       xx=[xx stst2.x];
	     end;
	     delta=delta-sys_deri(xx,hcli.parameter,0,[],[]);
	     for t=1:n_tau
	       delta=delta-sys_deri(xx,hcli.parameter,t,[],[])*exp(-lambda(i)*tau(t));
	     end;
	     delta=delta';
	[eigvec,eigval]=eig(delta);
	[minval pos]=min(abs(diag(eigval)));
	w(:,i)=eigvec(:,pos);
      end;
    end;
    
    hcli.v=v;
    hcli.w=w;
    
    hcli.alpha=hcli.v\(hcli.profile(:,1)-hcli.x1);
    hcli.epsilon=norm(hcli.alpha);
    hcli.alpha=hcli.alpha/hcli.epsilon;    
    
    return;
    
