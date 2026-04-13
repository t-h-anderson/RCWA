function [Z,Tn,R,E] ...
    = RCWA(inmk0, radtheta, epsf, nmdz, varargin)
arguments
    inmk0    (1,1) double {mustBeReal, mustBePositive}
    radtheta (1,1) double {mustBeReal, mustBeFinite}
    epsf     (:,:) double
    nmdz     (1,:) double {mustBeReal, mustBePositive}
    varargin
end
if size(epsf, 1) ~= length(nmdz)
    error('RCWA:sizeMismatch', ...
        'epsf must have one row per layer: expected %d rows, got %d.', ...
        length(nmdz), size(epsf, 1));
end
% [Z, Tn, R, e_x, e_y, e_z, h_x, h_y, h_z] = RCWA(epsmat, epsmat_recip, nmdz, (opt) varargin)
% Z is the Z matrix from EMSW
% Tn the tranfter matrix for each mode through each slice
% R the reflected wave
% E is the resulting electric field
% H is the resulting magnetic field

loc=DefaultLoc;
loc = ApplyVarargin(loc,varargin);
Ns = length(nmdz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start calculation
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Decompose light wave vectors (for incident, reflected and transmitted) - in air with refractive index nsa, grating
% with period nmLx
%radtheta = rad_from_deg(loc.degtheta);
inmk0x = loc.nsa*inmk0*sin(radtheta)+(-loc.Nt:loc.Nt)*2*pi/loc.nmLx;
inmk0z = sqrt(loc.nsa^2*inmk0^2-inmk0x.^2);

% Create Y matrices depending on polarization (ESW - eqns. 2.128 and 2.129)
if (loc.pol == 1) % nonzero right columns of Y matrices
    Ye_inc = - diag(inmk0z/(loc.nsa*inmk0));
    Ye_ref = - Ye_inc;
    
    Yh_inc = - diag(ones(2*loc.Nt+1,1));
    Yh_ref = Yh_inc;
    
else % nonzero left columns of Y matrices
    Ye_inc = loc.nsa*diag(ones(2*loc.Nt+1,1));
    Ye_ref = Ye_inc;
    
    Yh_inc = - diag(inmk0z/(loc.nsa*inmk0));
    Yh_ref = - Yh_inc;
end

Y_inc = [Ye_inc; Yh_inc]; % Y matrix for incident and transmitted field
Y_ref = [Ye_ref; Yh_ref];

Z  = cell(Ns+1,1);
Z{Ns+1} = Y_inc; % Eqn. 2.421

% Initialise matrices
Vn = cell(Ns+1,1);
inmGn = cell(Ns+1,1);
X_Upper = cell(Ns+1,1);
inmG_Upper = cell(Ns+1,1);

NullMat=zeros(2*loc.Nt+1);

inmKX=diag(inmk0x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     The implementation of the stable algorithm begins
% (The notation is taken from Electromagnetic Surface Waves -
%  T. Mackay, J. Polo, and A. Lakhtakia, p75-77)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% E=(Ex,Ey,Ez)=(Ex,0,Ez)+(0,Ey,0) and H=(Hx,Hy,Hz) = (Hx,0,Hz) +(0,Hy,0)
% Since E and H are orthogonal and by linearity, you can solve Maxwell
% for either E1=(Ex,0,Ez) and H1 = (0,Hy,0)
% or E2 = (0,Ey,0) and H2 = (Hx,0,Hz).
%(Ex,0,Ez)  and (0,Hy,0) corresponds to P-polaraization
%(0,Ey,0) and  (Hx,0,Hz) corresponds to S-polaraization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epst = zeros(Ns+1, 2*loc.Nt+1, 2*loc.Nt+1);
invepst = zeros(Ns+1, 2*loc.Nt+1, 2*loc.Nt+1);
epsz = 0;
for n = Ns:-1:1
    % Select the optical permittivity, and the inverse, for the slice
    
    if ~isequal(epsz,epsf(n,:)) || n==Ns
        epsz = epsf(n,:);
        % epsz_recip = epsmat_recip(:,n);
        
        if sum(epsz==0) >= 4*loc.Nt
            epst_temp = epsz(2*loc.Nt+1)*diag(2*loc.Nt+1);
            invepst_temp = (1/epsz(2*loc.Nt+1))*eye(2*loc.Nt+1);
                       
            epst(n,:,:) = epst_temp;
            invepst(n,:,:) = invepst_temp;
            
            if (loc.pol == 1)
                p_tr = sqrt(inmk0 - inmk0x.*inmk0x/(epsz(2*loc.Nt+1)*inmk0)); % p[n]
                p_bl = sqrt(inmk0 * epsz(2*loc.Nt+1));    % p
                
                pd_p = diag(p_tr/p_bl);
                p_pd = diag(p_bl./p_tr);
            else
                p_tr = sqrt(-inmk0);
                p_bl = sqrt(- inmk0*epsz(2*loc.Nt+1) + (1./inmk0)*inmk0x.*inmk0x);
                
                pd_p = diag(p_tr./p_bl);
                p_pd = diag(p_bl/p_tr);
            end
            
            V =  [pd_p -pd_p; eye(2*loc.Nt+1) eye(2*loc.Nt+1)];
            inmG = [p_bl*p_tr -p_bl*p_tr];
            invV = 0.5*[p_pd eye(2*loc.Nt+1) ; -p_pd eye(2*loc.Nt+1)];
            
            [V, order] = sortmat(V,inmG);
            invV=invV(order,:);
            inmG = inmG(order);
            
            
            
        else
            
            epst_temp= toeplitz(epsz(2*loc.Nt+1:-1:1),epsz(2*loc.Nt+1:4*loc.Nt+1));             % Eqn. 2.114
            invepst_temp = inv(epst_temp);
                       
            epst(n,:,:) = epst_temp;
            invepst(n,:,:) = invepst_temp;
               
            if (loc.pol == 1)
                inmP14 = inmk0*eye(2*loc.Nt + 1) - (1/inmk0)*inmKX*invepst_temp*inmKX;
                %inmP41 = inmk0*inv_epsf; % Why is this used?
                inmP41 = inmk0*epst_temp;
                
                inmP = [NullMat, inmP14; inmP41, NullMat];
            else
                inmP23 = - inmk0*eye(2*loc.Nt + 1);
                inmP32 = - inmk0*epst_temp + (1./inmk0)*inmKX*inmKX;
                
                inmP = ([NullMat, inmP23; inmP32, NullMat]);
            end
            
            % decompose P into form given in 2.138, order by the eigen values and
            
            [V, inmGD] = eig(inmP);
            inmG=diag(inmGD);
            [V, order] = sortmat(V,inmG);
            inmG = inmG(order);
            
            
        end
    else
         epst(n,:,:) = epst_temp;
         invepst(n,:,:) = invepst_temp;
    end
    
    
    if  sum(epsz==0) >= 4*loc.Nt
        X = invV*Z{n+1};
    else
        X = V\Z{n+1}; % X = inv(V).Z{n+1}
    end
    
    % store V and G in Vn and Nn for current layer
    Vn{n+1} = V;
    inmGn{n+1} = inmG;
    
    % Split G into upper and lower diagonal matrices
    inmG_Upper{n + 1} = inmG(1:2*loc.Nt + 1);
    inmG_Lower = inmG(2*loc.Nt + 2:4*loc.Nt + 2);
    
    inmexpG_U = diag(exp( 1i.*nmdz(n).*inmG_Upper{n + 1}));
    inmexpG_L = diag(exp(-1i.*nmdz(n).*inmG_Lower));
    
    % Can we make it analytic down to here?
    
    
    % Define the X matrix as in Eqn. 2.144
    X_Upper{n + 1} = X(1:2*loc.Nt + 1, :);
    
    % Define the U matrix as in Eqn. 2.146
    U = inmexpG_L*(X(2*loc.Nt + 2 : (4*loc.Nt + 2), :)*(X_Upper{n + 1}\inmexpG_U));
    
    % Calculate the next Z using Eqn. 2.145
    Z{n} = V*[eye(2*loc.Nt + 1); U];
end

% We can now calculate Tn{i}

% A vector describing the excited modes of the incident field.
A=zeros(2*loc.Nt+1,1);  % Y and Z are for the correct pol only, so problem dims are halved
A(loc.Nt+1) = 1;

% Eqn. 2.149 - [T0; R] = [Z(0), - Y_ref] * Y_inc * A

T0R = [Z{1}, - Y_ref]\(Y_inc*A); % Why is this the inverse...??

T0 = T0R(1 : 2*loc.Nt + 1, :);
R = T0R(2*loc.Nt + 2 : (4*loc.Nt + 2), :);

% Build Transmission Vectors
Tn=cell(Ns+1,1);
Tn{1} = T0;


for n = 1: Ns
    inmG_Upper = inmGn{n+1}(1: 2*loc.Nt + 1);
    Tn{n+1} = X_Upper{n+1}\(diag(exp(1i*nmdz(n)*inmG_Upper))*Tn{n}); % Eqn. 2.143
end


% Calculate the electric fields
% Create the f vector which contains [e_x; e_y; h_x; h_y]
f = zeros(Ns,4*loc.Nt+2);
e_x = zeros(Ns,2*loc.Nt+1);
e_y = zeros(Ns,2*loc.Nt+1);
e_z = zeros(Ns,2*loc.Nt+1);
%h_x = zeros(Ns,2*loc.Nt+1);
h_y = zeros(Ns,2*loc.Nt+1);
%h_z = zeros(Ns,2*loc.Nt+1);

E = zeros(3, Ns, 2*loc.Nt+1);
%H = zeros(3, Ns, 2*loc.Nt+1);

for n = 1: Ns
    
   % f(n,:) = 0.5*(Z{n}*Tn{n}+Z{n+1}*Tn{n+1}); % Eqn. 2.143
    f(n,:) = Z{n+1}*Tn{n+1};
    
    if(loc.pol == 1) % p-polarisation
        e_x(n,:) = f(n,1:2*loc.Nt+1);
        h_y(n,:) = f(n,2*loc.Nt+2:4*loc.Nt+2);
        %eps = toeplitz(epsf(2*loc.Nt+1:-1:1,n),epsf(2*loc.Nt+1:4*loc.Nt+1,n));
        %inv_eps = toeplitz(inv_epsmat(2*loc.Nt+1:-1:1,i),inv_epsmat(2*loc.Nt+1:4*loc.Nt+1,i));
        %e_z(i,:) = - 1./(im_from_inm(inmk0)*const.c) *  const.eps0 * inv(eps) * im_from_inm(inmKX)*(h_y(i,:)).';
        e_z(n,:) = -1/inmk0 *  squeeze(invepst(n,:,:)) * (inmKX)*(h_y(n,:)).';
        
        
        %        h_x(i,:) = 0*e_x(i);
        %        e_y(i,:) = 0*h_y(i);
        %        h_z(i,:) = 0*e_z(i);
    elseif(loc.pol == 0) % s-polarisation
        e_y(n,:) = f(n,1:2*loc.Nt+1);
        %        h_x(i,:) = f(i,2*loc.Nt+2:4*loc.Nt+2);
        % h_z(i,:) = (im_from_inm(inmk0)/(const.c*const.mu0)) * inmKX*(e_y(i,:)).';
        %        h_z(i,:) = im_from_inm(inmk0) * inmKX*(e_y(i,:)).';
        
        
        %        e_x(i,:) = 0*h_x(i);
        %        h_y(i,:) = 0*e_y(i);
        %        e_z(i,:) = 0*h_z(i);
    end
end

E(1,:,:) = e_x;
E(2,:,:) = e_y;
E(3,:,:) = e_z;

% H(1,:,:) = h_x;
% H(2,:,:) = h_y;
% H(3,:,:) = h_z;
