classdef Geometry
    % Vertical layering of the cell.
    %
    % Replaces BuildMaterial.m. From a Config, builds:
    %   nmz     - vector of z-sample positions (top of cell first; the
    %             original fliplr'd ordering is preserved)
    %   nmdz    - vector of slice thicknesses, aligned with nmz
    %   mat_cat - categorisation per z-slice:
    %               0 = mirror, 1 = grating, 2 = junction,
    %               3 = window, 4 = air
    %
    % Layer order from substrate (z = 0) upwards:
    %     mirror | grating | junction (p-i-n) | window | air

    properties (SetAccess = immutable)
        cfg     rcwa.Config
        nmz     (1,:) double
        nmdz    (1,:) double
        mat_cat (1,:) double
        nmgtop  (1,1) double    % z-coordinate of top of grating region
        nmLz    (1,1) double    % total junction thickness (nmLp+nmLi+nmLn)
    end

    methods
        function obj = Geometry(cfg)
            cfg = cfg.verify();
            obj.cfg = cfg;

            nmLz_ = cfg.nmLp + cfg.nmLi + cfg.nmLn;
            obj.nmLz = nmLz_;
            obj.nmgtop = cfg.nmda + cfg.nmLm;

            % Per-region z samples (cell midpoints, in ascending z order).
            [zm, dm] = obj.regionPoints(0, cfg.nmLm,  cfg.Nm, cfg.nmLm > 0);
            [zg, dg] = obj.regionPoints(cfg.nmLm, cfg.nmda, cfg.Ng, cfg.nmda > 0);
            zJ_start = cfg.nmLm + cfg.nmda;
            [zJ, dz] = obj.regionPoints(zJ_start, nmLz_, cfg.Nz, nmLz_ > 0);
            zW_start = zJ_start + nmLz_;
            [zw, dw] = obj.regionPoints(zW_start, cfg.nmLw, cfg.Nw, cfg.nmLw > 0);
            zA_start = zW_start + cfg.nmLw;
            [za, da] = obj.regionPoints(zA_start, cfg.nmLair, cfg.Nair, cfg.nmLair > 0);

            nmz_  = [zm, zg, zJ, zw, za];
            nmdz_ = [zeros(size(zm))+dm, zeros(size(zg))+dg, ...
                     zeros(size(zJ))+dz, zeros(size(zw))+dw, ...
                     zeros(size(za))+da];

            % Original flips so traversal in RCWA goes top->bottom.
            obj.nmz  = fliplr(nmz_);
            obj.nmdz = fliplr(nmdz_);

            obj.mat_cat = obj.categorise(obj.nmz, cfg);
        end
    end

    methods (Static, Access = private)
        function [zs, dz] = regionPoints(z0, L, N, has_region)
            if has_region && N >= 1
                dz = L / N;
                zs = linspace(z0 + 0.5*dz, z0 + L - 0.5*dz, N);
            else
                dz = 0;
                zs = [];
            end
        end

        function mat_cat = categorise(nmz, cfg)
            % Reproduce the original "halved sign trick" categorisation.
            nmLz = cfg.nmLp + cfg.nmLi + cfg.nmLn;

            mirror = sign(nmz - cfg.nmLm);
            mirror = 0.5 * (-mirror + abs(mirror));

            grating = sign(nmz - (cfg.nmLm + cfg.nmda)) + mirror;
            grating = 0.5 * (-grating + abs(grating));

            junction = sign(nmz - (nmLz + cfg.nmda + cfg.nmLm)) + mirror + grating;
            junction = 0.5 * (-junction + abs(junction));

            window = sign(nmz - (nmLz + cfg.nmda + cfg.nmLm + cfg.nmLw)) ...
                     + mirror + junction + grating;
            window = 0.5 * (-window + abs(window));

            air = sign(nmz - (nmLz + cfg.nmda + cfg.nmLm + cfg.nmLw + cfg.nmLair)) ...
                  + mirror + junction + grating + window;
            air = 0.5 * (-air + abs(air));

            mat_cat = 0*mirror + 1*grating + 2*junction + 3*window + 4*air;
        end
    end
end
