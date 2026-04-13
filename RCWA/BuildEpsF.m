function [epsf]=BuildEpsF(nmz, epsxz, loc)
arguments
    nmz   (1,:) double {mustBeReal, mustBeFinite}
    epsxz (:,:) double
    loc   (1,1) struct
end
if size(epsxz, 2) ~= length(nmz)
    error('BuildEpsF:sizeMismatch', ...
        'epsxz must have one column per z-slice: expected %d columns, got %d.', ...
        length(nmz), size(epsxz, 2));
end

%%%%%%%%%%%%%%%%%%%%%%%%Optimization%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

homo = zeros(length(nmz),1);

for i = 1:length(nmz)
    homo(i) = all(epsxz(:,i)==epsxz(1,i));
end

% if(loc.dmat==1) % Dielectric material in the grating
%     glasseps = (glass('nmlambda', nmlambda)+(10^-3)*1i).^2;
%     loc.epsd = glasseps;     % permittivity of the glass
%    
% end
% 
% if(loc.wmat==1) % Dielectric material in the window
%     glasseps = (glass('nmlambda', nmlambda)+(10^-3)*1i).^2;
%     loc.epsw = glasseps;     % permittivity of the glass
% end
% 
% if(loc.gmat == 1) % Metallic material in the grating
%     silvereps = (Silver('nmlambda', nmlambda).^2);
%     loc.epsm = silvereps;    % permittivity of the circle
% end

% Pass the material function, argument of position here, to the solver.
% Solved once nmz is populated in BuildEps
% if(loc.material == 1 && loc.Nz ~=0)
%     loc.epsJ = aSiHGC('eVEg', Eg(mat_cat==2).', 'nmlambda', nmlambda);
% elseif(loc.material ==2)
%     loc.epsJ = Faryad(nmz, nmlambda);
% end

% Build elist
% Begin building permittivity (eps) vectors here
Nz = size(nmz);
epsf=zeros(Nz(2), 4*loc.Nt+1);

% epsf(2*loc.Nt+1, mat_cat == 4) = loc.nsa;
% 
% epsf(2*loc.Nt+1, mat_cat == 3) = loc.epsw;
% 
% epsf(2*loc.Nt+1, mat_cat == 2) = loc.epsJ;
% 
% epsf(2*loc.Nt+1, mat_cat == 0) = loc.epsm;

epsf( homo == 1, 2*loc.Nt+1) = epsxz(1,homo==1);

if (loc.ftype == 1)
    epsf(homo ~= 1, :) = epsn(loc, nmz(homo ~= 1));
else
    fftresult = gepszfft(nmz(homo ~= 1), loc);
    epsf(homo ~= 1, :) = fftresult;
end
