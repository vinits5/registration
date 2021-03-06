clearvars -except a;
clc;

% Parameters.
sensor_file = 'sensor_50pts_normalized_noisy_data.csv';
model_file = 'model_100pts_normalized.csv';
sigma = 0.01;
rot_range = 30;     % In degrees about ground truth rotation. (euler angles)
trans_range = 0.2;  % In percentage about ground truth rotation.
show_datas = false; % Display sensor & model data in 3d plots.
no_rand_transforms = 10000; % No of random transformations to be used to evaluate the obj function.
itr_nos = 40;       % No of optimization iterations.
list_Omin = [];     % Store the min objective function after each iteration.
Omin = 0; tmin = [0;0;0]; euler_min = [0;0;0];  % Just initialization of variables.
range_reduction_factor = 0.9;   % Reduce the range to select random tranformations

helper = HelperRegistration();   % Create Object of Helper Class.

% Read Sensor (Nx3) and Model (Mx3) Data.
[sensor,model] = helper.read_data(sensor_file,model_file);

% Ground Truth Values.
gt_rotation = [80;20;0];        % Euler Angles (in degrees)
gt_translation = [0;0.1;0.1];   % Ground Truth translation
[ref_rotation,ref_translation] = helper.sensor2model_Transform(gt_rotation,gt_translation);  % Find Rotation (3x1) & Translation (3x1) which can convert sensor data to model data.
%ref_rotation = gt_rotation;
%ref_translation = gt_translation;
dlmwrite('results.txt',ref_rotation','-append','delimiter',' ','roffset',0);
dlmwrite('results.txt',ref_translation','-append','delimiter',' ','roffset',0);
disp(ref_rotation');
disp(ref_translation');

sensor = helper.apply_transformation(sensor,gt_rotation,gt_translation);             % Apply rotation & translation on model data to get sensor data.
gt_Oval = helper.call_obj_emicp(sensor,model,ref_rotation,ref_translation,sigma);      % Minimum value of Obj function at ground truth transfromation.
disp('test');
disp(gt_Oval);
dlmwrite('results.txt',gt_Oval,'-append','delimiter',' ','roffset',0);
% Display Sensor & Model Data.
if show_datas==true
    scatter3(sensor(1:end,1),sensor(1:end,2),sensor(1:end,3));
    figure();
    scatter3(model(1:end,1),model(1:end,2),model(1:end,3));
end

% Used for seeding random number generators.
%a = rng;
%save('seeder','a');
%load('seeder.mat');
%rng(a);

tic;    % To calculate time required for the computation.
for idx = 1:itr_nos
    fprintf('idx no: %d\n',idx);
    fprintf('No. Random pts: %d\n',no_rand_transforms);
    % Call optimization function.
    [Omin, euler_min, tmin, check_prev_cond] = helper.optimization(no_rand_transforms,ref_rotation,ref_translation,rot_range,trans_range,sigma,sensor,model,Omin,idx);
    % Display Results in command window and store in .txt file.
    helper.display_results(idx,Omin,euler_min,tmin,check_prev_cond);
    
    % Update your reference rotation & translation with the rotation & translation given by optimizer.
    ref_translation = tmin;
    ref_rotation = euler_min;
    
    % Update the random search region about referance transformation.
    if rem(idx,1)==0
        % Update the rotation range about the ref_rotation in Euler Angles.
        rot_range = rot_range*range_reduction_factor;
        % Update the translation range about the ref_translation.
        trans_range = trans_range*range_reduction_factor;
        %rot_range = max(rot_range,3);
        %trans_range = max(trans_range,0.01);
        disp(rot_range);
        disp(trans_range);
    end
    list_Omin = [list_Omin,Omin];
end
elasped_time=toc;       % To calculate time required for computation.
plot(list_Omin);        % Plot the objective function. 
disp(elasped_time);