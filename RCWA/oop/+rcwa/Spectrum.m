classdef Spectrum
    % AM1.5G solar spectrum tabulator (W m^-2 nm^-1).
    %
    % Loads AM15G.csv once on construction and provides interpolation at
    % an arbitrary set of wavelengths. Loader uses readmatrix when
    % available (R2019a+), otherwise falls back to csvread.

    properties (SetAccess = immutable)
        nmlambda_table (:,1) double
        flux_table     (:,1) double
    end

    methods
        function obj = Spectrum(csvPath)
            if nargin < 1 || isempty(csvPath)
                csvPath = 'AM15G.csv';
            end
            if exist('readmatrix', 'file') == 2
                data = readmatrix(csvPath);
            else
                data = csvread(csvPath);
            end
            obj.nmlambda_table = data(:, 1);
            obj.flux_table     = data(:, 2);
        end

        function W = spectralIrradiance(obj, nmlambda)
            % Linear interpolation of the table at the given wavelengths.
            W = interp1(obj.nmlambda_table, obj.flux_table, nmlambda);
        end
    end
end
