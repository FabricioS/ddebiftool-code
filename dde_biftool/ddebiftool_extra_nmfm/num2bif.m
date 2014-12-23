function biftype = num2bif(index)
% Convert an index to a bifurcation type
% Return number of supported types on 'count'

switch index
   case 0
      biftype = 'stst';
   case 1
      biftype = 'hopf';
   case 2
      biftype = 'fold';
   case 3
      biftype = 'psol';
   case 4
      biftype = 'hcli';
   case 5
      biftype = 'genh';
   case 6
      biftype = 'hoho';
   case 7
      biftype = 'zeho';
   case 'count'
      biftype = 7;
   otherwise
      error('NUM2BIF: unknown bifurcation index = %g', index);
end


end

