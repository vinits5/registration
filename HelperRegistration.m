% Helper class to define useful functions for registration.

classdef HelperRegistration
    properties
        objective_function_number = 0;      % Assign the function number to use that objective function.
        sigma_x = 0.01^2; sigma_y = 0.02^2; sigma_z = 0.015^2;     % Diagonal Elements of covariance matrix for noise.
    end
    methods
        % Generate sensor data (randomly select points from model data).
        function generate_sensor_data(~,file_name,num_sensor_pts,display_data)
            % Arguments-->
            % file_name:        name of model file                          (ex. 'model_100pts.csv')
            % num_sensor_pts:   No of points in sensor dataset.             (No. of pts x 3)
            % display_data:     plots sensor & model data                   boolean
            
            model = csvread(file_name);
            % sensor = [];                                                   % Array to store sensor points.
            shuffled_indices = randperm(100);                                % Shuffled indices to get the sensor data.
            sensor = zeros(num_sensor_pts,3);
            for i = 1:num_sensor_pts
                %sensor = [sensor; model(randi(size(model,1)),1:end)];       % Stores random points from model data.
                sensor(i,1:end) = model(shuffled_indices(i),1:end);
            end
            dlmwrite('correspondences.csv',shuffled_indices,'delimiter',',');   % Store Correspondences in a file.
            csvwrite('sensor_50pts_normalized.csv',sensor);                 % Store sensor data to .csv file.
            if display_data
                scatter3(model(1:end,1),model(1:end,2),model(1:end,3));     % Plot model data.
                figure();
                scatter3(sensor(1:end,1),sensor(1:end,2),sensor(1:end,3));  % Plot sensor data.
            end
        end
        
        % Read Sensor and Model files.
        function [sensor,model] = read_data(~,sensor_file,model_file)
            % Arguments-->
            % sensor_file                   Name of file containing sensor points in '.csv' format. (String)
            % model_file                    Name of file containing model points in '.csv' format. (String)
            
            % Ouput-->
            % sensor:                       Sensor Points; float (Nx3)
            % model:                        Model Points; float (Mx3)
            
            sensor = csvread(sensor_file);
            model = csvread(model_file);
        end
        
        % Find Value of Objective Function (EMICP).
        % objective_function_number = 0
        % Find Value of Objective Function (EMICP).
        function value = call_obj_emicp(obj,sensor_pts,model_pts,euler,trans,sigma)
            % Arguments-->
            % sensor_pts:   Matrix of Sensor Points                             Nx3
            % model_pts:    Matrix of Model Points                              Mx3
            % euler:        Euler Angles                                        3x1
            % trans:        Translation Matrix                                  3x1
            % sigma:        Sigma Value for Normal Distribution (Variance)      1x1
            
            % Output-->
            % value         Evaluated value of objective function               1x1
            value = 0;
            Rot = obj.euler2rotation(euler);                                        % Convert euler angles to rotation matrix.
            for i = 1:size(sensor_pts,1)                                            % Loop for Sensor pts.
                temp_i = 0;
                for k = 1:size(model_pts,1)                                           % Loop for Model pts.
                    temp_vec = (model_pts(k,1:end)'-Rot*sensor_pts(i,1:end)'-trans);  % si-R*mk-t := temp_vec
                    temp_k = temp_vec'*temp_vec;                                      % matrix_transpose*matrix := temp_k
                    temp_k = exp(temp_k*(-1/(2*sigma*sigma)));                        % exp((-1/2*sgima_square)*temp_k)
                    temp_i = temp_i + temp_k;                                         % Sum of temp_k for all model points
                end
                const = size(model_pts,1)*sqrt(2*pi*sigma*sigma);                     % 1/(Nm*sqrt(2*pi*sigma_square))
                temp_i = log(temp_i/const);
                if isinf(temp_i)
                    temp_i=-300;                                                    % ITS A FLOATING POINT ISSUE. WILL DISCUSS LATER.
                end
                value = value + temp_i;                                             % Sum of all temp_i for sensor points.
            end
            value = -value;
        end
    
        % Objective Function 1
        % objective_function_number = 1
        function value = call_icp_obj1(obj,sensor_pts,model_pts,euler,trans)
            % Arguments-->
            % sensor_pts    Matrix of Sensor Points             Nx3
            % model_pts     Matrix of Model Points              Mx3
            % euler         Euler Angles (in degrees)           3x1
            % trans         Translation Matrix                  3x1
            
            % Output-->
            % value         Evaluated value of objective function   1x1
            
            Rot = obj.euler2rotation(euler);                    % Convert euler angles to rotation matrix.
            temp_mat = model_pts'-Rot*sensor_pts';              % temp_mat = mi-R*si    (3xN)
            temp_mat = temp_mat-trans;                          % temp_mat = temp_mat-t = mi-R*si-t (3xN)
            temp_val = temp_mat'*temp_mat;                      % temp_val (NxN)
            value = sum(diag(temp_val));                        % sum of diagonal elements.
        end
    
        % Objective Function 2
        % objective_function_number = 2
        function value = call_icp_obj2(obj,sensor_pts,model_pts,euler,trans)
            % Arguments-->
            % sensor_pts    Matrix of Sensor Points             Nx3
            % model_pts     Matrix of Model Points              Mx3
            % euler         Euler Angles (in degrees)           3x1
            % trans         Translation Matrix                  3x1
            
            % Output-->
            % value         Evaluated value of objective function   1x1
            
            Rot = obj.euler2rotation(euler);                                    % Convert euler angles to rotation matrix.
            temp_mat = model_pts'-Rot*sensor_pts';              % temp_mat = mi-R*si    (3xN)
            temp_mat = temp_mat-trans;                          % temp_mat = temp_mat-t = mi-R*si-t (3xN)
            
            Sigma = diag([obj.sigma_x,obj.sigma_y,obj.sigma_z]);    % Define covariance matrix
            B = Rot*Sigma*Rot';                                     
            B = inv(B);                                             % Bi = inv(R'*Cov/Sigma*R)
            temp_val = temp_mat'*B*temp_mat;                        % (mi-R*si-t)'*Bi*(mi-R*si-t)    (NxN)
            value = sum(diag(temp_val));                            % sum of diagonal elements.
        end
            
        % Convert euler angles to Rotation matrix.
        function Rot = euler2rotation(~,euler)
            % Arguments-->
            % euler(1): Rotation about x-axis in degrees.
            % euler(2): Rotation about y-axis in degrees.
            % euler(3): Rotation about z-axis in degrees.
            
            % Output-->
            % Rot:      Rotation Matrix with 'xyz' as euler angles.
            Rot_x = [1, 0, 0; 0, cosd(euler(1)), -sind(euler(1)); 0, sind(euler(1)), cosd(euler(1))];
            Rot_y = [cosd(euler(2)), 0, sind(euler(2)); 0, 1, 0; -sind(euler(2)), 0, cosd(euler(2))];
            Rot_z = [cosd(euler(3)), -sind(euler(3)), 0; sind(euler(3)), cosd(euler(3)), 0; 0, 0, 1];
            Rot = Rot_x*Rot_y*Rot_z;
        end
    
        % Find rotation and translation which can be applied to sensor pts to get model pts.
        function [gt_euler,gt_translation] = sensor2model_Transform(~,euler,translation)
            % Arguments-->
            % euler            Rotation Angles (in degrees) to convert from model pts to sensor pts. (3x1)
            % translation      Translation to convert from model pts to sensor pts. (3x1)
            
            % Outputs-->
            % gt_euler         Rotation Angles (in degrees) to convert from sensor pts to model pts. (3x1)
            % gt_translation   Translation to convert from sensor pts to model pts. (3x1)
            
            Rot_euler = eul2rotm(euler'*(pi/180),'XYZ');                        % Rotation matrix from euler angles.
            Transformation = zeros(4,4);
            Transformation(4,4) = 1;
            Transformation(1:3,1:3) = Rot_euler;
            Transformation(1:3,4)= translation;                                 % Model to Sensor Transformation Matrix.
            Transformation = inv(Transformation);                               % Sensor to Model Transformation Matrix.
            temp_euler = rotm2eul(Transformation(1:3,1:3),'XYZ')*(180/pi);      % Euler Angles (in degrees) (Sensor to Model) (1x3)
            gt_euler = temp_euler';                                             % Euler Angles (in degrees) (Sensor to Model) (3x1)
            gt_translation = Transformation(1:3,4);                            % Translation (Sensor to Model) (3x1)
        end
            
        % Generate random rotation and translation about a particular euler angle and translation.
        function [euler_rotation_angles,translation] = random_transformation(~,ref_rotation,ref_translation,rot_range,trans_range)
            % Arguments-->
            % ref_rotation:     Reference Rotation (in Degrees) (Euler Angles) (Find random rotation about these values)
            % ref_translation:  Reference Translation (Find random translation about these values)
            % rot_range:        Range of Rotation (in Degrees) (Euler Angles)
            % trans_range:      Range of Translation (%)
            
            % Output-->
            % euler_rotation_angles:    Random Euler Angles in given range around reference rotation angles.
            % translation:              Random translation in given range around reference translation.
            euler_rotation_angles = 2*rot_range*rand(3,1)-rot_range+ref_rotation;           % [ref_rotation-rot_range, ref_rotation+rot_range]
            
            % [ref_translation-trans_range*ref_translation, ref_translation+trans_range*ref_translation]
            translation = 2*trans_range*ref_translation.*rand(3,1)-trans_range*ref_translation+ref_translation;
            
            % Check if the ref_translation is zero. 
            for i = 1:size(translation,1)
                if translation(i)==0
                    translation(i) = 2*trans_range*rand()-trans_range;  % [-trans_range, trans_range]
                end
            end
        end
        
        % Apply Ground Truth Transformation to Sensor Data.
        function transformed_sensor_data = apply_transformation(obj,sensor_data,gt_rotation,trans)
            % Arguments-->
            % sensor_data:      Randomly sampled data from model. (Nx3)
            % gt_rotation:      Apply this rotation on model data to generate sensor data. (euler angle 3x1)
            % trans:            Translate model data to generate sensor data. (3x1)
            
            % Output-->
            % transformed_sensor_data       Sensor Point cloud rotated and translated as per ground truth value. (Nx3)
            
            Rot = obj.euler2rotation(gt_rotation);                          % Euler angles to rot matrix
            transformed_sensor_data = Rot*sensor_data';                     % Apply rotation(3x3) to sensor data(3xN).
            transformed_sensor_data = transformed_sensor_data + trans;      % Apply translation(3x1) to sensor data(3xN).
            transformed_sensor_data = transformed_sensor_data';             % Transpose of sensor data. (Nx3)
        end
        
        % Optimization Algorithm
        function [Omin, euler_min, tmin, check_prev_cond] = optimization(obj,no_random_pts,ref_rotation,ref_translation,rot_range,trans_range,sigma,sensor,model,Omin_prevItr,idx)
            % Arguments-->
            % no_random_pts                     No of random points to be chosen.
            % ref_rotation                      Euler angles obtained in previous iteration (3x1)
            % ref_translation                   Translation obtained in previous iteration (3x1)
            % rot_range                         Range about ref_rotation to choose random euler angles.
            % rot_trans                         Range about ref_translation to choose random translation values.
            % sigma                             Sigma value for Normal Probability Distribution to calculate obj value.
            % sensor                            Sensor Points (Nx3)
            % model                             Model Points (Mx3)
            % Omin_prevItr                      Minimum Objective Value in last iteration.
            % idx                               Iteration number which is going on.
            
            % Outputs-->
            % Omin                              Minimum Objective Value in current iteration.
            % euler_min                         Euler angles which has given minimum objective value (Omin).
            % tmin                              Translation vector which has given minimum objective value (Omin).
            % check_prev_cond                   True: If Omin of previous iteration is less than Omin of current iteration else, False.
            
            all_obj_val = zeros(1,no_random_pts);       % Store obj values for all random transformations in given range. (1xRandom_Pts)
%             Rmin = zeros(3,loop_nos);
            all_euler_min = [];                     % Store all Random Euler Angles generated 
            all_tmin = [];                          % Store all Random Translations generated.
%             tmin = zeros(3,loop_nos);
            parfor i = 1:no_random_pts              % parfor to evaluate obj function in parallel.
                % Generate random euler_rotation_angles (3x1) and translation (3x1) in a given range to evaluate obj function.
                [euler_rotation_angles,translation] = obj.random_transformation(ref_rotation,ref_translation,rot_range,trans_range);
                if obj.objective_function_number == 0
                    all_obj_val(i) = obj.call_obj_emicp(sensor,model,euler_rotation_angles,translation,sigma);    % Find value of emicp obj function.
                elseif obj.objective_function_number == 1
                    all_obj_val(i) = obj.call_icp_obj1(sensor,model,euler_rotation_angles,translation);           % Find value of obj function 1
                elseif obj.objective_function_number == 2
                    all_obj_val(i) = obj.call_icp_obj2(sensor,model,euler_rotation_angles,translation);           % Find value of obj function 2
                end
                all_euler_min = [all_euler_min, euler_rotation_angles];                         % Store euler rotation angles. (generated randomly)
                all_tmin = [all_tmin, translation];                                             % Store randomly generated translation.
                
                % Only useful if for loop is used. (Saves the memory)
%                 if i == 1
%                     Omin = val; Rmin = rotation; tmin = translation;
%                 else
%                     if val<Omin
%                         Omin = val; Rmin = rotation; tmin = translation;
%                     end
%                 end
            end
            [obj_val,idx_Omin] = min(all_obj_val);                                  % Find minimum objective value.
            
            if (obj_val>Omin_prevItr && idx ~= 1)                                   % Compare the min obj value with obj value obtained in previous iteration.
                % Return all the results obtained in previous iteration.
                Omin = Omin_prevItr;
                euler_min = ref_rotation;
                tmin = ref_translation;
                check_prev_cond = true;
            else
                % Return the new min obj value and rotation and translation associated with it.
                Omin = obj_val;
                euler_min = all_euler_min(1:end,idx_Omin);
                tmin = all_tmin(1:end,idx_Omin);
                check_prev_cond = false;
            end
        end
        
        % Display Results in the command window and store in results.txt file.
        function display_results(~,idx,Omin,euler_min,tmin,check_prev_cond)
            % Arugments-->
            % idx               Iteration number
            % Omin              Minimum objective value after iteration number idx.
            % euler_min         Euler rotation angles (in degrees) to result Omin.
            % tmin              Translation which results Omin.
            % check_prev_cond   True: if the Omin belongs to previous iteration. Flase: if Omin belongs to current iteration.
            
            fprintf('Min Obj value: %f\n',Omin);
            fprintf('Rotation: ');
            disp(euler_min);
            fprintf('Translation: ');
            disp(tmin);
            fprintf('Is Prev Obj Min Less: ');
            disp(check_prev_cond);
            
            % Store these values in results.txt file.
            if idx == 1
                dlmwrite('results.txt',idx,'-append','delimiter',' ','roffset',0);
            else
                dlmwrite('results.txt',idx,'-append','delimiter',' ','roffset',1);
            end
            dlmwrite('results.txt',Omin,'-append','delimiter',' ','roffset',0);
            dlmwrite('results.txt',euler_min','-append','delimiter',' ','roffset',0);
            dlmwrite('results.txt',tmin','-append','delimiter',' ','roffset',0);
        end
    
        % Add Guassian Noise to the sensor data.
        function sensor_noise_pts = add_noise(~,sensor_pts)
            % Argument-->
            % sensor_pts        Points in Sensor point cloud (Nx3)
            
            % Output-->
            % sensor_noise_pts  Points with Guassian Noise added in Sensor Point Cloud. (Nx3)
            
            sensor_noise_pts = zeros(size(sensor_pts));
            for i=1:size(sensor_pts,1)
                mu = [0,0,0];
                sigma_x2 = 0.01^2;        % sigma(x^2)
                sigma_y2 = 0.02^2;        % sigma(y^2)
                sigma_z2 = 0.015^2;        % sigma(z^2)
                % sigma = [sigma_x2,0,0;0,sigma_y2,0;0,0,sigma_z2];
                sigma = diag([sigma_x2, sigma_y2, sigma_z2]);           % Covariance Matrix.
                rand_noise = mvnrnd(mu,sigma);                                   % Function to sample random guassian noise.
                sensor_noise_pts(i,:) = sensor_pts(i,:)+rand_noise;              % Add noise to sensor point.
            end
        end
        
        % Find table using obj functions and their obtained minimum parameters.
        function table = find_obj_table(obj,sensor_pts,model_pts,corr_model_pts,no_obj_func,euler,trans,sigma)
            % Argument-->
            % sensor_pts                Sensor data (Nsx3)
            % model_pts                 Model data (Nmx3)
            % corr_model_pts   Model data (which corresponds to sensor data) (Nsx3)
            % no_obj_func               No of objective functions(N)         (scalar)
            % euler                     Euler Angles for each obj function which will produce min Obj value (in degree) (3xN).
            % trans                     Translation for each obj function which will produce min Obj value (3xN)
            % sigma                     Sigma for the covariance matrix.
            
            % Output-->
            % table             A matrix with all the obj values for each rotation and translation for each function (NxN)
            
            if no_obj_func==size(euler,1)                       % Just to ensure that there is a set of euler angles corresponding to each objective function.
                table = zeros(no_obj_func,no_obj_func);         % Create an empty matrix to store all objective values
                for i=1:size(euler,1)
                    table_col = zeros(no_obj_func,1);
                    % Call a function which will evaluate all obj functions for given orientation & translation.
                    for j = 0:no_obj_func-1
                        table_col(j+1,1) = obj.evaluate_obj_value(j,sensor_pts,model_pts,corr_model_pts,euler(:,i),trans(:,i),sigma);
                    end
                    table(:,i) = table_col;
                end
            else
                disp('Dimension Error: No of objective functions and Set of Euler angles are not equal');
            end
        end
        
        % Evaluate objective function to help find_obj_table function to find a column of the table using given rotation & translation.
        function obj_value = evaluate_obj_value(obj,func_no,sensor_pts,model_pts,corr_model_pts,euler,trans,sigma)
            % Arguments-->
            % func_no               Define numbers for each function.
            % sensor_pts                Sensor data (Nsx3)
            % model_pts                 Model data (Nmx3)
            % corr_model_pts   Model data (which corresponds to sensor data) (Nsx3)
            % euler                     Euler Angles for each obj function which will produce min Obj value (in degree) (3x1).
            % trans                     Translation for each obj function which will produce min Obj value (3x1)
            % sigma                     Sigma for the covariance matrix.
            
            % Output-->
            % Objective Value of given function for given euler angles and translations.
            
            if func_no == 0
                obj_value = obj.call_obj_emicp(sensor_pts,model_pts,euler,trans,sigma);
            elseif func_no == 1
                obj_value = obj.call_icp_obj1(sensor_pts,corr_model_pts,euler,trans);
            else
                obj_value = obj.call_icp_obj2(sensor_pts,corr_model_pts,euler,trans);
            end
        end
        
    end
end