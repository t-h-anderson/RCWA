function epsz = epsn(loc, nmz)
arguments
    loc (1,1) struct
    nmz (1,:) double {mustBeReal, mustBeFinite}
end
% Only allow gratings that do not contain semiconductor material
if(loc.nmda < loc.nmLg)
    loc.nmda = loc.nmLg;
end

% Build grating eps
n = -2*loc.Nt:1:2*loc.Nt;
nz = size(nmz);
epsz = zeros(nz(2), 4*loc.Nt+1);


mloc = 0.5*(1 + sign(loc.nmLm - nmz)); % Array has value 1 in mirror range
gloc = 0.5*(1 + sign(loc.nmLg + loc.nmLm - nmz)) - mloc; % Array has value 1 in mirror and metallic grating region
dloc = -mloc-gloc+1;

sz=size(epsz(gloc == 1, n~=0)); % Size of grating region
one = ones(sz(1),1);

if (loc.type==0)   % Square Grating  
      % Zero mode is the average permittivity
      epsz(mloc == 1, n==0) = loc.epsm;
      epsz(gloc == 1, n==0) = (loc.epsd*(1-loc.zeta) + loc.epsm*loc.zeta);
      epsz(dloc == 1, n==0) = loc.epsd;
      
      % Fourier modes occur only in the grating region
      epsz(gloc == 1, n~=0) = one*((loc.epsm-loc.epsd)*(sin(n(n~=0)*pi.*loc.zeta))./(pi*n(n~=0)));     
      epsz(gloc == 0, n~=0) = 0;
elseif (loc.type == 1)
      % Zero mode is the average permittivity
      
      width = 2*acos((nmz(gloc==1)-loc.nmLm)/loc.nmLg) * loc.zeta/pi;
          
      epsz(mloc == 1, n==0) = loc.epsm;
      epsz(gloc == 1, n==0) = (loc.epsd*(1-width) + loc.epsm*width);
      epsz(dloc == 1, n==0) = loc.epsd;
      
      % Fourier modes occur only in the grating region
      %size(n(n~=0))
      %size(width)
      %size(epsz(n~=0, gloc == 1) )
      %size(kron(n(n~=0).',width))
      epsz(gloc == 1, n~=0) = ((loc.epsm-loc.epsd)*(sin(kron(n(n~=0), pi*width.')))./(kron(n(n~=0),pi*one)));     
      epsz(gloc == 0, n~=0) = 0;
elseif (loc.type == 2)
      % Zero mode is the average permittivity
      
      width = (1-((nmz(gloc==1)-loc.nmLm)/loc.nmLg)) * loc.zeta;
          
      epsz(mloc == 1, n==0) = loc.epsm;
      epsz(gloc == 1, n==0) = (loc.epsd*(1-width) + loc.epsm*width);
      epsz(dloc == 1, n==0) = loc.epsd;
      
      % Fourier modes occur only in the grating region
      %size(n(n~=0))
      %size(width)
      %size(epsz(n~=0, gloc == 1) )
      %size(kron(n(n~=0).',width))
      epsz(gloc == 1, n~=0) = ((loc.epsm-loc.epsd)*(sin(kron(n(n~=0), pi*width.')))./(kron(n(n~=0),pi*one)));     
      epsz(gloc == 0, n~=0) = 0;
else
    error('Exact Grating Fourier Not Implemented. Add method in epsn.m.');
end

%epsz(2:2:end,:) = conj(epsz(2:2:end,:));

