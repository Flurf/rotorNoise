classdef Propeller
    properties
        airfoil %airfoil type, NACA rules
        propfile
        A   % area of disk rotor
        R   % total radius of prop
        R0  %radius at root from rotation oriin
        n   % total number of blades
        beta0   %root pitching angle
        gamma   %twist ratio
        r   %array of radius
        c   %array of chords
        beta %array of angles of pitch
        

    end
    methods
        
        function obj = Propeller(R,R0,n,rn,chord,AoAroot,twist,airfoiltype)
            if nargin == 8
                obj.R = R;
                obj.R0 = R0;
                obj.n = n;
                obj.r = linspace(obj.R0,obj.R,rn);           %constant radius spacing
                obj.c = zeros(1,length(obj.r));                %constant chord
                obj.c = obj.c + chord;
                obj.beta0 = AoAroot;
                obj.gamma = twist;
                obj.airfoil = airfoiltype;

                obj.A = obj.rotorArea(obj.R,obj.R0);
    
                obj.beta = obj.twist(obj.r,obj.R,obj.R0,obj.beta0,obj.gamma);
            end
        end

        function f = readProp(filename)
            %TODO
        end
        function t = twist(obj,r,R,R0,beta0,gamma)
            t = zeros(1,length(r));
            for i=1:length(r)
                t(i) = beta0 - gamma*(r(i) - R0)/(R - R0);
            end
        end
        function a = rotorArea(obj,R,R0)
            a = pi*(R)^2 - pi*(R0)^2;
        end
    end
end