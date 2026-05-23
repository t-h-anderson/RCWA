classdef Constants
    % Physical constants in SI units.
    properties (Constant)
        mu0  = 4*pi*1e-7;            % Permeability of free space [H/m]
        eps0 = 8.854e-12;            % Permittivity of free space [F/m]
        eta0 = 376.73031346177;      % Impedance of free space   [Ohm]
        c    = 3e8;                  % Speed of light            [m/s]
        h    = 6.6260696e-34;        % Planck constant           [J s]
        q    = 1.6021766208e-19;     % Elementary charge         [C]
    end
end
