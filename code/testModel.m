function testModel(skeleton_num,muscle_tendon_parameters_num,casadiFun,vizualizationFun)
Muscle_activation = .15 ; 
%% test of kinematics 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % skeleton spatial configuration 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q1 = 0; % x
q2 = 0; % y
q3 = 0; % z
q4 = 0; % alpha hip
q5 = 0; % alpha knee
q6 = -25; % alpha ankle

% conversion deg to rad
q5 = (q5/180)*pi; q6 = (q6/180)*pi; 

q_num = [q1,q2,q3,q4,q5,q6];

    % plot skeleton 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Origin_num,Insertion_num,~,Markers_num] = casadiFun.ForwardKinematics([q_num,skeleton_num]);

plotmodel(Origin_num, Insertion_num, Markers_num)

    % mtu length
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mtu_length_num  = full(casadiFun.getUMTLength([q_num,skeleton_num]));

tibialis_length_num = mtu_length_num(1);
soleus_length_num = mtu_length_num(2);
gastrocnemius_length_num = mtu_length_num(3);

disp(['tibials muscle-tendon unit length: ', num2str(tibialis_length_num),' m'])
disp(['soleus muscle-tendon unit length: ', num2str(soleus_length_num),' m'])
disp(['gastrocnemius muscle-tendon unit length: ', num2str(gastrocnemius_length_num),' m'])

    % mtu moment arm
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mtu_moment_arm_num  = full(casadiFun.getMomentArm([q_num,skeleton_num]));

tibialis_moment_arm_num = mtu_moment_arm_num(1,end);
soleus_moment_arm_num = mtu_moment_arm_num(2,end);
gastrocnemius_moment_arm_num = mtu_moment_arm_num(3,end);

disp(['tibials muscle-tendon unit moment arm: ', num2str(tibialis_moment_arm_num),' m'])
disp(['soleus muscle-tendon unit moment arm: ', num2str(soleus_moment_arm_num),' m'])
disp(['gastrocnemius muscle-tendon unit moment arm: ', num2str(gastrocnemius_moment_arm_num),' m'])

    %% test of muscle dynamics 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% exemple of gastrocnemius 
Optimal_fiber_length_2 = muscle_tendon_parameters_num(3);
Pennation_angle_at_muscle_optimal_fiber_length_2 = muscle_tendon_parameters_num(6);
Maximal_isometric_muscle_force_2 = muscle_tendon_parameters_num(9);
Tendon_slack_length_2 = muscle_tendon_parameters_num(12);
muscle_tendon_parameters_2 = muscle_tendon_parameters_num([3,6,9,12]);

% test tendon 
% tendon_current_length = Tendon_slack_length_2 + (Tendon_slack_length_2*.03);
% plotTendonForcelength(vizualizationFun,Tendon_slack_length_2,Maximal_isometric_muscle_force_2,tendon_current_length)

% test muscle 
% Fiber_length = Optimal_fiber_length_2 + (Optimal_fiber_length_2*.3);
% plotMuscleForcelength(vizualizationFun,Optimal_fiber_length_2,Maximal_isometric_muscle_force_2,Muscle_activation,Fiber_length)

    %% test rooted variables (signle muscle) - exemple soleus 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
known_num = [Muscle_activation,gastrocnemius_length_num , muscle_tendon_parameters_2];
% unknown  = vertcat(tendonLength(1), fiberLength(1), pennationAngle(1)) ;
% as x0 (xstart)
unknown  = [Tendon_slack_length_2, Optimal_fiber_length_2, Pennation_angle_at_muscle_optimal_fiber_length_2] ;

rooted_variables= full(casadiFun.equilibrateMuscleTendonSingleMuscle(unknown,known_num));

tendon_length_rooted = rooted_variables(1);
fiber_length_rooted = rooted_variables(2);
pennation_angle_rooted = rooted_variables(3);
pennation_deg = (pennation_angle_rooted/pi)*180;


%vizualization
plotTendonForcelength(vizualizationFun,Tendon_slack_length_2,Maximal_isometric_muscle_force_2,tendon_length_rooted)
plotMuscleForcelength(vizualizationFun,Optimal_fiber_length_2,Maximal_isometric_muscle_force_2,Muscle_activation,fiber_length_rooted)

disp(['rooted fiber length: ', num2str(fiber_length_rooted),' m'])
disp(['rooted tendon length: ', num2str(tendon_length_rooted),' m'])
disp(['rooted pennation angle: ', num2str(pennation_deg),' deg'])


muscleActiveForce = full(Muscle_activation* vizualizationFun.getNormalizedMuscleActiveForce(fiber_length_rooted,Optimal_fiber_length_2)*Maximal_isometric_muscle_force_2);
musclePassiveForce = full(vizualizationFun.getNormalizedMusclePassiveForce(fiber_length_rooted,Optimal_fiber_length_2)*Maximal_isometric_muscle_force_2);
muscleTotalForce = muscleActiveForce + musclePassiveForce;

disp(['muscle passive force: ', num2str(musclePassiveForce),' N'])
disp(['muscle active force: ', num2str(muscleActiveForce),' N'])
disp(['muscle total force: ', num2str(muscleTotalForce),' N'])

% error 
residuals = full(casadiFun.equilibriumErrorSingleMuscle(unknown',known_num)); % residuals
disp(['residuals g5: ', num2str(residuals(1))])
disp(['residuals g6: ', num2str(residuals(2))])
disp(['residuals g7: ', num2str(residuals(3))])

end 