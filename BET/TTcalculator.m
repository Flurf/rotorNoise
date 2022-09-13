classdef TTcalculator
    properties
        density = 1.225;
        viscosity = 0.0000148;
        c0 = 340;       %reference sound speed
        Tmargin = 0.1;  %margin of T - Texp, decimal value


        Propeller
        T
        Texp        %Thrust expected, used to extimate induced velocity
        RPM
        psi        %array of azimuth, in radians
        V           %total speed
        Vp          %speed on rotation plane
        Vf          %forward flight speed
        Vind        %induced speed
        Re          %array of effective reynolds numbers (induced + rotational speed + forward speed)
        AoA
        M
        cl          %array of lift coeffs
        cd          %array of drag coeffs
        L
        Tl
        cp           
    end


    methods
        function obj = TTcalculator(Propeller,RPM,Vf,deltapsi,Texp,Vind)
            if nargin
                obj.Propeller = Propeller;
                obj.Texp = Texp;

                obj.psi = linspace(0,360,360/deltapsi);
                obj.RPM = RPM;
                obj.Vf = Vf;
                if ~exist('Vind','var')
                     obj.Vind = obj.inducedVelocity(obj.Propeller.r);
                else
                    % i can give induce velocity from file. It has to be of
                    % the same length of Prop.r and obv. sampled in the
                    % same points as r.
                    obj.Vind = Vind;
                end
                obj.Vp = obj.planeVelocity(obj.Vf, obj.RPM, obj.Propeller.r,obj.Propeller.n);
                obj.V = obj.velocity(obj.Vp, obj.Vind,obj.Propeller.n);
                obj.AoA = obj.trueAoa(obj.Vp,obj.Vind,obj.Propeller.n);
                obj.Re = obj.reynolds(obj.Propeller.c,obj.V,obj.Propeller.n);
                obj.M = obj.mach(obj.V,obj.Propeller.n);
                [obj.cl,obj.cd,obj.L,obj.Tl, obj.cp] = obj.lift(obj.Propeller.airfoil,obj.Re,obj.AoA,obj.M,obj.Propeller.n);
                obj.T = obj.thrust(obj.Propeller.r,obj.psi,obj.Propeller.n,obj.Tl);
            end
        end

        function V = velocity(obj,Vp,Vind,n)
           V = zeros(1,length(Vind),length(obj.psi),n);
           for k = 1:n
               for i=1:length(obj.psi)
                   for j = 1:length(Vind)
                        V(1,j,i,k) = sqrt((Vind(1,j))^2 + (Vp(1,j,i,k))^2);
                        if Vp(1,j,i,k)    < 0
                            V(1,j,i,k) = 0;
                            %if velocity enters from Trailing edge don't
                            %consider that section, no lift produced 
                        end       
                   end
               end
           end
        end


        function Vp = planeVelocity(obj,Vf, RPM,r,n)
            Vrot = zeros(1,length(r));  %rotational speed
            Vplane = zeros(1,length(r),length(obj.psi),n);  %rotation plane speed (forward + rotational) for each blade n
            
            for i = 1:length(Vrot)
                Vrot(i) = (pi*RPM/30)*r(i);
            end

            for j = 1:length(obj.psi)
                for i = 1:length(Vplane(1,:,1,1))  
                    for k = 1:n
                        Vplane(1,i,j,k) = Vrot(1,i) - Vf*sin( (obj.psi(j) + (k-1)*(360/n))*(pi/180) );           
                    end
                end
            end

            Vp = Vplane;


        end


        function vind = inducedVelocity(obj,r)
                vind = zeros(1,length(r));
                vind = vind + sqrt(obj.Texp/(obj.Propeller.A*2*obj.density));
                %TODO other complex induced speed estimations are possible 

        end

        function t = trueAoa(obj,Vp,Vind,n)
            AoA = zeros(1,length(obj.Propeller.beta),length(obj.psi),n);

            for k = 1:n
                for i = 1:length(obj.psi)
                    for j = 1:length(obj.Propeller.beta)
                        if obj.V(1,j,i,k) ~= 0
                            AoA(1,j,i,k) = obj.Propeller.beta(1,j) - atan(Vind(1,j)/Vp(1,j,i,k))*(360/(2*pi));
                        else
                            AoA(1,j,i,k) = 0;
                        end
                        if AoA(1,j,i,k) < -15
                            AoA(1,j,i,k) = 0;
                        end
                    end
                end
            end
            t = AoA;
        end

        function r = reynolds(obj,c,V,n)
            r = zeros(1,length(c),length(obj.psi),n);
            for k = 1:n
                for i = 1:length(obj.psi)
                    for j = 1:length(c)
                        r(1,j,i,k) = obj.density*c(1,j)*V(1,j,i,k)/obj.viscosity;
                    end
                end
            end
        end

        function M = mach(obj,V,n)
            M = zeros(1,length(obj.Propeller.r),length(obj.psi),n);
            for k = 1:n
                for i = 1:length(obj.psi)
                    for j = 1:length(obj.Propeller.r)
                        M(1,j,i,k) = V(1,j,i,k)/obj.c0;
                    end
                end
            end
        end
        
        function [cl,cd,L,Tl,cpr] = lift(obj,airfoil, Re, AoA, M,n)

            cl = zeros(1,length(Re),length(obj.psi));
            cd = zeros(1,length(Re),length(obj.psi));
            cpr = zeros(1,length(Re),length(obj.psi));
            Ltop = 0;  %peak of load after wich parabolic decrease
            Tltop = 0;
            Rtop = 0;
            for k = 1:n
                for i = 1:length(obj.psi)
                    for j = 1:length(obj.Propeller.r)
                        
                        %[P,F] = xfoil(airfoil,AoA(:,j,i),0,0,'oper/iter 800', 'panels n 100')
                        if (obj.Propeller.r(j) < 1.5)
                            [P,F] = xfoil(airfoil,AoA(:,j,i,k),0,0,'oper/iter 800', 'panels n 100');
                        else
                            [P,F] = xfoil(airfoil,AoA(:,j,i,k),0,M(:,j,i,k),'oper/iter 800', 'panels n 100');
                        end
                        %[P,F] = xfoil(airfoil,AoA(:,j,i),Re(:,j,i),M(:,j,i),'oper/iter 800', 'panels n 100')
                        cl(1,j,i,k) = P.CL;
                        cd(1,j,i,k) = P.CD;
                    
                         xd= F.xcp(1:length(F.xcp)/2);
                         xv= F.xcp(length(F.xcp)/2 + 1:length(F.xcp));
                        
                         cpd= F.cp(1:length(F.cp)/2);
                         cpv= F.cp(length(F.cp)/2 + 1:length(F.cp));
                         cp = abs(trapz(xd,cpd)) - abs(trapz(xv,cpv));
                         cpr(1,j,i,k) = cp;
                         L(1,j,i,k) = 0.5 * obj.density *  obj.Propeller.c(j) * (obj.V(1,j,i,k))^2 * cp ;
                         Tl(1,j,i,k) = 0.5 * obj.density *  obj.Propeller.c(j) *(obj.V(1,j,i,k))^2 * (P.CL* cos(AoA(1,j,i,k) * (2*pi/360)) -P.CD* sin(AoA(1,j,i,k) * (2*pi/360)));
                         %parabolic decrease of loads for r>0.95r
                         if obj.Propeller.r(j) > 0.97 * obj.Propeller.R
                            if Ltop == 0 
                                Ltop = L(1,j,i,k);
                                Tltop = L(1,j,i,k);
                                Rtop = obj.Propeller.r(j);
                            end

                            L(1,j,i,k) = Ltop * (1-((obj.Propeller.r(j) - Rtop)/(obj.Propeller.R - Rtop))^2);
                            Tl(1,j,i,k) = Tltop * (1-((obj.Propeller.r(j) - Rtop)/(obj.Propeller.R - Rtop))^2);
                         end
                    end
                end
            end
        end
        
        function T = thrust(obj,r,psi,n,Tl)
            T = zeros(1,length(psi));
            for i = 1:length(psi)
                for k = 1:n
                                    
                    T(1,i) = T(1,i) + trapz(r, Tl(1,:,i,k));
                end
            end
        end
    end
 end