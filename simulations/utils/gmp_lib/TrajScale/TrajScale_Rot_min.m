%% TrajScale_Rot_min class
%
% For generalizing a 3-DoF trajectory to new target/initial positions using 
% a rotation and scaling base on
% scale type 1 from 
% 'A novel DMP formulation for global and frame independent spatial scaling in the task space'
% DOI: 10.1109/RO-MAN47096.2020.9223500
% 
% scaling matrix: 
% Ks = ( ||g - y0|| / ||gd - yd0|| ) * {rotation matrix that aligns gd - yd0 with g - y0 using the minimum angle of rotation}
% 
% where g, y0 are the new target and initial positions
%      gd, yd0 are the demonstrated target and initial positions
%

classdef TrajScale_Rot_min < TrajScale
    
    methods (Access = public)
        
        %% Constructor.
        %  @param[in] n_dof: degrees of freedom.
        function this = TrajScale_Rot_min()
            
            this@TrajScale(3);
            
        end

    end
    
    methods (Access = public) % Abstract implementations

        function scale_type = getScaleType(this)

            scale_type = TrajScale.ROT_MIN_SCALE;

        end

    end

    methods (Static)
        function R = custom_axang2rotm(axang)
            % Custom implementation of axang2rotm
            % axang = [axis_x, axis_y, axis_z, angle]
            k = axang(1:3) / norm(axang(1:3));  % Normalize axis
            theta = axang(4);
            K = [0, -k(3), k(2); k(3), 0, -k(1); -k(2), k(1), 0];  % Skew-symmetric matrix
            R = eye(3) + sin(theta)*K + (1 - cos(theta))*K^2;
        end
    end
    
    methods (Access = protected)
       
        % ------------------------------------------
        
        function sc = calcScaling(this, Y0, Yg)
            
            this.Y0 = Y0;
            this.Yg = Yg;
            
            nd = this.Ygd - this.Y0d;  nd = nd/norm(nd);
            n = this.Yg - this.Y0;  n = n/norm(n);
            dot_n_nd = dot(n,nd);
            if (abs(abs(dot_n_nd) - 1) < 1e-14)
                R = eye(3,3);
            else
                k = cross(nd,n);
                theta = acos(dot_n_nd);
                R = TrajScale_Rot_min.custom_axang2rotm([k' theta]);
            end
            
            sc = R*norm(this.Yg - this.Y0)/norm(this.Ygd - this.Y0d);
            
        end
        
        function sc = calcInvScaling(this)
        
            nd = this.Ygd - this.Y0d;  nd = nd/norm(nd);
            n = this.Yg - this.Y0;  n = n/norm(n);
            dot_n_nd = dot(n,nd);
            if (abs(abs(dot_n_nd) - 1) < 1e-14)
                R = eye(3,3);
            else
                k = cross(nd,n);
                theta = acos(dot_n_nd);
                R = TrajScale_Rot_min.custom_axang2rotm([k' theta]);
            end
            sc = R'*norm(this.Ygd - this.Y0d)/norm(this.Yg - this.Y0);
            
        end

        % ------------------------------------------
        
    end
    
end
