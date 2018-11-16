classdef MCML_model2 < handle
    %   This code constructs light transport in multi-layers model
    %   Coded by Zheyuan
    
    properties
        iterations = 10000; % iteration times
        ni=1; % reflection index of input
        nt=1.37; % reflection index of output
        slice = 200; % the size of total dividence
        threshold_dead = 0.00001; % the dead threshold for photon package
        prob_survive = 0.01; % the possibility of survive
        x=0;   % position of photon package along x
        y=0;   % position of photon package along y
        z=0;   % position of photon package along z(depth direction)
        vx=0;  % direction of photon package along x
        vy=0;  % direction of photon package along y
        vz=1;  % direction of photon package along z(depth direction)
        weight=1; % weight of each photon package
        out_bool=0;  % whether the photon package is outside
        dead_bool=0; % the survive condition of photon package
        distance=0;
        hit_bool=0;
        
        rs=0;
        final_result; % to record the final weight absorption matrix
        temp_result;
        R_final_result; % to record the final weight reflectance
        T_final_result; % to record the final weight transmit
        Z_final_result;
        total_result;
        
        ua; % absorption coefficient
        us; % scatter coefficent
        ut; % interaction coefficient
        g;  % anisotropy <cos(theta)>
        d;  % the depth of tissue
        r;  % size of tissue:-r:r
        db=0; % the hit distance
    end
    
    methods
        function obj = MCML_model2(ua, us, g, d, r, ni, nt, iterations, slice)
            %   MCML_model: Construct an instance of this class
            %   Basic parameter of MCML
            obj.ua = ua;
            obj.us = us;
            obj.g = g;
            obj.d = d;
            obj.ni = ni;
            obj.nt = nt;
            obj.r = r;
            obj.iterations = iterations;
            obj.slice = slice;
            % obj.ut = ua + (1 - g)*us;
            obj.ut = ua + us;
            
            obj.total_result = zeros(d*slice,r*slice, r*slice);
            obj.temp_result = zeros(d*slice,r*slice, r*slice);
            obj.final_result = zeros(d*slice,r*slice); % to record the final weight absorption matrix
            obj.R_final_result = zeros(r*slice,1); % to record the final weight reflectance
            obj.T_final_result = zeros(r*slice,1); % to record the final weight transmit
            obj.Z_final_result = zeros(d*slice,1); % to record the Z direction
        end
        
        function obj=initial_package(obj)
            obj.x=0;   % position of photon package along x
            obj.y=0;   % position of photon package along y
            obj.z=0;   % position of photon package along z(depth direction)
            obj.vx=0;  % direction of photon package along x
            obj.vy=0;  % direction of photon package along y
            obj.vz=1;  % direction of photon package along z(depth direction)
            obj.weight=1-((obj.ni-obj.nt)/(obj.ni+obj.nt))^2; % weight of each photon package
            obj.out_bool=0;  % whether the photon package is outside
            obj.dead_bool=0; % the survive condition of photon package
            obj.distance=0;
        end 
        
        function obj=scatter(obj)
            %scatter calculate the scattered vector for the photon
            %   For updating 
            theta_sample = rand;
            
            if obj.g==0
                cos_theta = 2*theta_sample - 1;
            else
                cos_theta = 0.5*(1+obj.g^2- ((1-obj.g^2)/(1-obj.g+2*obj.g*theta_sample))^2)/obj.g;
                
                if abs(cos_theta)>1
                    cos_theta = sign(cos_theta);
                end
                
            end
            
            phi_sample = rand;
            phi = 2*pi*phi_sample;
            
            sin_theta = sqrt(1-cos_theta^2);
            cos_phi = cos(phi);
            sin_phi = sin(phi);
            
            if abs(obj.vz)<0.9999
                c = sqrt(1-obj.vz^2);
                vx2 = sin_theta*(obj.vx*obj.vz*cos_phi - obj.vy*sin_phi) / c + obj.vx*cos_theta;
                vy2 = sin_theta*(obj.vy*obj.vz*cos_phi + obj.vx*sin_phi) / c + obj.vy*cos_theta;
                vz2 = - c*sin_theta*cos_phi + obj.vz*cos_theta;
            else
                vx2 = sin_theta*cos_phi;
                vy2 = sin_theta*sin_phi;
                vz2 = sign(obj.vz)*cos_theta;
            end
            
            obj.vx = vx2;
            obj.vy = vy2;
            obj.vz = vz2;
            
        end
        
        function obj=move(obj,distance)
            obj.x = obj.x + obj.vx*distance;
            obj.y = obj.y + obj.vy*distance;
            obj.z = obj.z + obj.vz*distance;
        end
        
        function obj=is_hit_boundary(obj)
            if obj.vz < 0
                obj.db = -obj.z/obj.vz;
%             elseif obj.vz > 0
%                 obj.db = (obj.d-obj.z)/obj.vz;
            else
                obj.db = Inf;
            end
            
            obj.hit_bool = (obj.db*obj.ut <= obj.distance);
            
        end
        function obj=hit_boundary(obj)
           
            ai=acos(abs(obj.vz));
            if obj.ni<obj.nt 
                angle_full_reflect = asin(obj.ni/obj.nt);
                if ai<=angle_full_reflect
                    at=asin(sin(ai)*obj.nt/obj.ni);
                    Rai=1/2*((sin(ai-at)/sin(ai+at))^2 + ...
                        (tan(ai-at)/tan(ai+at))^2);
                else
                    Rai=1;
                end
            else
                at=asin(sin(ai)*obj.nt/obj.ni);
                Rai=1/2*((sin(ai-at)/sin(ai+at))^2 + ...
                    (tan(ai-at)/tan(ai+at))^2);
            end
            
            random_sample = rand;
            if random_sample <= Rai
                 obj.vz=-obj.vz;
            else
                 obj.vx=obj.nt*obj.vx/obj.ni;
                 obj.vy=obj.nt*obj.vy/obj.ni;
                 obj.vz=sign(obj.vz)*cos(at);
                 
                 obj.out_bool=1;
                 obj.dead_bool=1;
            end
            
        end
        
        function obj=absorption_record(obj, delta_weight)
            
            dis_x = abs(obj.x)*obj.slice;
            x_id = round(dis_x)*sign(obj.x)+1+obj.slice/2;
            % w1 = dis_x - floor(dis_x);
            
            dis_y = abs(obj.y)*obj.slice;
            y_id = round(dis_y)*sign(obj.y)+1+obj.slice/2;
            
            z_id = floor(obj.z*obj.slice)+1;
            
            if z_id>=1&&z_id<=(obj.d*obj.slice)
                obj.Z_final_result(z_id)=obj.Z_final_result(z_id) + delta_weight;
            end
            
            if z_id<=obj.d*obj.slice&&z_id>=1&&x_id<=obj.r*obj.slice&&x_id>=1&&y_id<=obj.r*obj.slice&&y_id>=1
                obj.total_result(z_id, x_id, y_id)=obj.total_result(z_id, x_id, y_id) + ...
                delta_weight;
            end
            
        end
        
        function obj=check_dead(obj)
            if obj.dead_bool~=1
                if obj.weight < obj.threshold_dead
                    obj.dead_bool = 1;
                    
                    random_sample = rand;
                    if random_sample < obj.prob_survive
                        obj.weight = obj.weight/obj.prob_survive;
                        obj.dead_bool = 0;
                    else
                        delta_weight=obj.weight;
                        obj.absorption_record(delta_weight);
                        obj.weight = 0;
                    end
                end
            end
        end
        function obj=simulate_once(obj)
            while obj.dead_bool==0
                
                if obj.distance==0
                    distance_sample = rand;
                    obj.distance = -log(distance_sample);
                end
                
                obj=obj.is_hit_boundary();
                
                if obj.hit_bool==1
                    obj=obj.move(obj.db);
                    obj.distance = obj.distance - obj.db*obj.ut;
                    obj=obj.hit_boundary();
                    
                    
                else
                    obj=obj.move(obj.distance/obj.ut);
                    
%                     if abs(obj.x)>=(obj.r/2)||abs(obj.y)>=(obj.r/2)||obj.z>obj.d
%                         obj.dead_bool=1;
%                         break;
%                     else
%                         delta_weight = obj.ua/obj.ut*obj.weight;
%                         obj.weight = obj.weight - delta_weight;
%                         obj=obj.absorption_record(delta_weight);
%                     end
                    delta_weight = obj.ua/obj.ut*obj.weight;
                    obj.weight = obj.weight - delta_weight;
                    obj=obj.absorption_record(delta_weight);
                    
                    obj.distance = 0;
                    obj=obj.scatter();
                    
                end
                
                obj=obj.check_dead();
            end
        end
        function obj=simulate(obj)
            rng('shuffle');
            for i=1:obj.iterations
                 obj=obj.initial_package();
                obj=obj.simulate_once();
            end
        end
        function imshow_absorption(obj)
            figure(1)
            z_absorption = obj.Z_final_result*obj.slice/(obj.iterations*obj.ua);
            semilogy((1:obj.slice)/obj.slice, z_absorption(1:obj.slice));
            title('Absorption Along Z Direction')
            xlabel('Z cm')
            ylabel('Absorption')
            figure(2);
            imagesc(log(squeeze(sum(obj.total_result,1))))
            title('Absorption in the surface')
            
            colormap jet;
            figure(3);
            imagesc(log(squeeze(sum(obj.total_result,2))))
            title('Absorption in the slice')
            
            colormap jet;
            s_log=log(z_absorption(round(obj.slice/4):obj.slice/2));
            x_log=(round(obj.slice/4):obj.slice/2)./obj.slice;
            p = polyfit(x_log', s_log, 1);
            sprintf('%f, %f',p(1), p(2))
            file_name=[num2str(obj.iterations),'_',num2str(obj.nt*100),'model.mat'];
            save(file_name,'obj');
        end
    end
end

