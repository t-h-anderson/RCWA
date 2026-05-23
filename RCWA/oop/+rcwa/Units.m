classdef Units
    % Unit conversions used throughout the RCWA pipeline.

    methods (Static)
        function m = m_from_nm(nm)
            m = nm * 1e-9;
        end

        function nm = nm_from_m(m)
            nm = m * 1e9;
        end

        function inv_m = im_from_inm(inv_nm)
            % Convert inverse nanometres to inverse metres: (1/nm) -> (1/m).
            inv_m = inv_nm * 1e9;
        end

        function rad = rad_from_deg(deg)
            rad = deg * pi/180;
        end

        function deg = deg_from_rad(rad)
            deg = rad * 180/pi;
        end

        function eV = eV_from_nm(nm)
            % Photon energy E = hc/lambda. Constant 1239.842 nm·eV is hc/q.
            eV = 1239.842 ./ nm;
        end

        function nm = nm_from_eV(eV)
            nm = 1239.842 ./ eV;
        end
    end
end
