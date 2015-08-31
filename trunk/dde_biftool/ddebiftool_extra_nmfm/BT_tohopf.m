function hopf=BT_tohopf(funcs,btpoint,freqs) %#ok<INUSD,INUSL>
hopf = struct('kind','hopf','parameter',btpoint.parameter,...
    'x',btpoint.x,'v',btpoint.q0,'omega',0);
if isfield(btpoint,'stability')
    hopf.stability=btpoint.stability;
end
end 
  