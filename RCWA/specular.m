function [rp0,rp1,tp0,tp1]=specular...
    (R,T, inmk0, radtheta, loc, varargin)
arguments
    R        (:,1) double
    T        (:,1) double
    inmk0    (1,1) double {mustBeReal, mustBePositive}
    radtheta (1,1) double {mustBeReal, mustBeFinite}
    loc      (1,1) struct
    varargin
end
if length(R) ~= length(T)
    error('specular:sizeMismatch', ...
        'R and T must have the same length: |R|=%d, |T|=%d.', ...
        length(R), length(T));
end

%%%%%%%%%%%%%%%%%%%%%%%%Optimization%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

loc = ApplyVarargin(loc,varargin);

% computations of transmission \
% coefficients, TP[i] contains the tranmission coefficients at the ith
% interface so TP[Ns+1] contains the coefficients of whole structure as
% their are Ns+1 interfaces

inmkxn = loc.nsa.*inmk0.*sin(radtheta)+(-loc.Nt:loc.Nt)*2*pi/loc.nmLx;
inmkzn = sqrt(loc.nsa^2.*inmk0^2-inmkxn.^2);

RP2 = abs(R).^2;    % Reflection amplitudes
TP2 = abs(T).^2;    % Transmission amplitudes


% The following block computes a diagonal matrix with real parts of \
% kz^(n) as its elements**)
kzr1 = real(inmkzn)/(loc.nsa*inmk0*cos(radtheta));
kzr2 = real(inmkzn)/(loc.nsa*inmk0*cos(radtheta));
RKZ1 = diag(kzr1);
RKZ2 = diag(kzr2);

%  (******************************* p incident, p reflected and \
% transmitted Pp ***************************************************

RPp = (RP2'*RKZ1).';
TPp = (TP2'*RKZ2).';
rp0 = RPp(loc.Nt+1);
rp1 = sum(RPp)-RPp(loc.Nt+1);
tp0 = TPp(loc.Nt+1);
tp1 = sum(TPp)-TPp(loc.Nt+1);
% (***************************************************************