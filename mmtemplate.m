function tL = mmtemplate(type,N)
 %
 % tL = mmtemplate(type,N)
 %
 %  Computes theoretical SQF responses for prototypical objects.
 %
 %     N - filter orders (from mmsteerablemap)
 %  type - one of: 'line', 'endpoint'
 %
 %
 % See also: mmsteerablemap
 
 % (c) Radim Sara (sara@cmp.felk.cvut.cz) FEE CTU Prague, 05 Jun 13
 
 switch lower(type)
  
  case 'line'
   if any(isodd(N))
    warning([mfilename, ':Inefficient'], 'Odd filter orders are not necessary for type "%s"',type)
   end
   
   tL = line_tmpl_n([0,pi], N);
   tL = tL./sqrt(tL*tL' - tL(1)^2);
 
  case 'endpoint'
   tL = line_tmpl_n(0, N);
 
  case 'edge'
   tL = wedge_tmpl_n([-pi/2,pi/2], N);

  case 'checker'
   d = pi/2;
   tL = wedge_tmpl_n([-d/2, d/2, pi-d/2, pi+d/2], N);
   %tL = checker_tmpl_n(N);
   % rotate by pi/4
   angle = pi/4;
   r = exp(1i*N*angle);
   tL = r.*tL;
  otherwise
   error([mfilename,':Unknown'], 'Unknown type %s',type)
 end
end

function tL = checker_tmpl_n(N)
 tL = sdiv(1i*(1+exp(1i*N*pi)).*sin(N*pi/4), N, 0);
end

function tL = wedge_tmpl_n(psi,N)
 tL = 0;
 for k=1:length(psi)/2
  c = exp(1i*N*psi(2*k-1)) - exp(1i*N*psi(2*k));
  tL = tL + c;   % line
 end
 tL = sdiv(tL,N,0);
end

function tL = line_tmpl_n(psi,N)
 tL = 0;
 for k=1:length(psi)
  tL = tL + exp(1i*N*psi(k));   % line
 end
end


