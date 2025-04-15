function [casadiFun,vizualizationFun,unknown_parameters,definition] = DeGrooteFunction()
% import Casadi 
[~, name] = system('hostname');
if strcmp(name(1:9), '942-27984')
    addpath('C:\Users\Stage\Desktop\Doctorat\Manip_Neuromusculoskeletal_Modeling\Casadi')
    elseif strcmp(name(1:8), '151302-1')
    addpath('C:\Users\amariani\Desktop\Thèse\Manip_Neuromusculoskeletal_Modeling\Casadi')
elseif strcmp(name(1:27), 'MacBook-Air-de-mickaelbegon')
    addpath('/Users/mickaelbegon/Downloads/casadi-3.6.3-osx64-matlab2018b/')
end
import casadi.*

%% Model definition 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % variables
    %%%%%%%%%%%
% "a" is neuromuscular activation (between 0 and 1) [3]
% "q" is the spatial skeleton configuration (x, y, z, theta_hip,theta_knee, theta_ankle) [6]
% "skeleton" is subject segment length and muscle insertion position. (segment_geometry, muscle_insertion, Local_ViaPoint_tibialis_anterior) [33]
% "muscleTendonParameters" : musculoskeletal parameters that are generally assumed not to change: maximal isometric muscle force (Fom), optimal fiber length (ℓom), tendon slack length (ℓst), and pennation angle at optimal fiber length (φo).

    % rooted
    %%%%%%%%%%%
% fiberLength [3]
% tendonLength [3]
% rootedvariables = [fiberLength, tendonLength] [6]

    % concatenated
    %%%%%%%%%%%
% "musculoskeletal_states" = [q, known_parameter] [39]
% "neuromusculoskeletal_state" = [a, q, known_parameter] [42]
% "all_states" = [neuromusculoskeletal_state, rooted_variales] [48]

                       %% 1. Musculoskeletal Geometry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 1.1 Topology 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Definition 
    %%%%%%%%%%%%%
% B : Bones
% J : Joins 
% df : degree of freedom
     
% Total df = 6 df
% --> J : Hip (4 df : Rotx= 0 ; Roty = 0 ; Rotz = theta_hip ; X = x ; Y = y; Z = z )
%   --> B : Thigh
%    --> J : knee (1 df : Rotx= 0 ; Roty = 0 ; Rotz = theta_knee ; X = 0 ; Y = 0 ; Z = 0 )
%       --> B : Leg 
%           --> J : Ankle (1 dF: Rotx= 0 ; Roty = 0 ; Rotz = theta_ankle ; X = 0 ; Y = 0 ; Z = 0 )
%                --> B : Talus 
%                   --> J : Subtalar (0 dF: Rotx= 0 ; Roty = 0 ; Rotz = 0 ; X = 0 ; Y = 0 ; Z = 0 )
%                        --> B : Calc 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nMuscles = 3 ; 
% NameMuscles = {'TibialisAnterior','Soleus','Gastrocnemius'} ; 
% nBones = 4 ; NameBones = {'Thigh','Leg','Talus','Calcaneus'} ; 
% nJoints = 4 ; NameJoints = {'Hip','Knee','Ankle','Subtalar'} ; 

    %% 1.2 Equation of the musculoskeletal model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.2.1 States (q's)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = SX.sym('x');
y = SX.sym('y');
z = SX.sym('z');
theta_hip = SX.sym('theta_hip');
theta_knee = SX.sym('theta_knee');
theta_ankle = SX.sym('theta_ankle');

q = vertcat(x,y,z,theta_hip, theta_knee, theta_ankle ); % skeleton configuration

% 1.2.2 Model parameters that are personalized by scaling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thigh = SX.sym('thigh',3);
leg = SX.sym('leg',3);
talus = SX.sym('talus',3);
clac = SX.sym('foot',3);

segment_geometry = vertcat(thigh, leg, talus ,clac);

% 1.2.3 Rototanslation matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R_0_thigh = Rototranslation_Rz(x,y,z, theta_hip);
R_thigh_leg = Rototranslation_Rz(thigh(1),thigh(2),thigh(3), -theta_knee);
R_leg_talus = Rototranslation_Rz(leg(1),leg(2),leg(3), -theta_ankle );
R_talus_calc = Rototranslation_Rz(talus(1),talus(2),talus(3), 0 );

% in R0
%R_0_leg = R_0_thigh * R_thigh_leg; 
%R_0_talus = R_0_leg * R_leg_talus; 
% R_0_calc = R_0_talus * R_talus_calc; 

% 1.2.4 Joint center
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Muscle insertions in local frame
Local_Insertion_tibialis_anterior = SX.sym('Local_Insertion_tibialis_anterior',3);
Local_Origin_tibialis_anterior = SX.sym('Local_Origin_tibialis_anterior',3);
Local_ViaPoint_tibialis_anterior = SX.sym('Local_ViaPoint_tibialis_anterior',3);

Local_Insertion_soleus = SX.sym('Local_Insertion_soleus',3);
Local_Origin_soleus = SX.sym('Local_Origin_soleus',3);
Local_Insertion_gastrocnemius = SX.sym('Local_Insertion_gastrocnemius',3);
Local_Origin_gastrocnemius = SX.sym('Local_Origin_gastrocnemius',3);

muscle_insertion = vertcat( ...
    Local_Origin_tibialis_anterior, ...
    Local_Origin_soleus, ...
    Local_Origin_gastrocnemius, ...
    Local_Insertion_tibialis_anterior, ...
    Local_Insertion_soleus, ...
    Local_Insertion_gastrocnemius ); %TODO: add viapoint

skeleton = vertcat(segment_geometry, muscle_insertion, Local_ViaPoint_tibialis_anterior);

% 1.2.5 Muscle origins and insertion in thigh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % rototranslation matrix 
R_thigh_talus = R_thigh_leg * R_leg_talus; 
R_thigh_calc = R_thigh_talus * R_talus_calc; 

    % joint center 
HJC = [0; 0; 0]; % Hip
KJC = R_thigh_leg(1:3, 4); % Knee
AJC = R_thigh_talus(1:3, 4); % Ankle
TJC = Rototranslate(R_thigh_calc, clac) ; % toe
CALC = R_thigh_calc(1:3, 4);  %calcaneum

    % insertion origine muscle 
Origin_tibialis_anterior = Rototranslate(R_thigh_leg, Local_Origin_tibialis_anterior);  
Insertion_tibialis_anterior = Rototranslate(R_thigh_calc, Local_Insertion_tibialis_anterior);
ViaPoint_tibialis_anterior = Rototranslate(R_thigh_leg, Local_ViaPoint_tibialis_anterior);

Origin_soleus = Rototranslate(R_thigh_leg,  Local_Origin_soleus) ;
Insertion_soleus = Rototranslate(R_thigh_calc, Local_Insertion_soleus);

Origin_gastrocnemius = Local_Origin_gastrocnemius;
Insertion_gastrocnemius = Rototranslate(R_thigh_calc, Local_Insertion_gastrocnemius) ;

Origin = horzcat(Origin_tibialis_anterior, Origin_soleus, Origin_gastrocnemius);
Insertion = horzcat(Insertion_tibialis_anterior, Insertion_soleus, Insertion_gastrocnemius);
ViaPoint = horzcat(ViaPoint_tibialis_anterior) ; 
Markers = horzcat(HJC, KJC, AJC, TJC,CALC); 

    %% 1.3 Casadi functions about model geometry/kinematics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
musculoskeletal_states = vertcat(q, skeleton) ; % all parameters of musculoskeletal model 

% ForwardKinematics
ForwardKinematics = Function('ForwardKinematics', ...
    {musculoskeletal_states}, {Origin, Insertion, ViaPoint, Markers}, ...
    {'musculoskeletal_states'}, {'Origin', 'Insertion', 'ViaPoint', 'Markers'}) ;

% umtLength : function to get UMT length of each muscle (TA, SOL, GAST)
% according to  the musculoskeletal model and his articular configutation 
umtLength = vertcat( sqrt(sum((Insertion_tibialis_anterior - ViaPoint_tibialis_anterior).^2)) + ...
                     sqrt(sum((ViaPoint_tibialis_anterior - Origin_tibialis_anterior).^2)) , ... % Tibialis
                     sqrt(sum((Insertion_soleus - Origin_soleus).^2)) , ... % Soleus
                     sqrt(sum((Insertion_gastrocnemius - Origin_gastrocnemius).^2)) ) ; % Gast

getUMTLength = Function('UMTLength', ...
    {musculoskeletal_states}, {umtLength}, ...
    {'musculoskeletal_states'}, {'umt_length'});

% momentArm : function to get Moment Arm of each muscle (TA, SOL, GAST)
% according to  the musculoskeletal model and his articular configutation 
momentArm = jacobian(umtLength, q);

getMomentArm = Function('MomentArm', ...
    {musculoskeletal_states}, {momentArm}, ...
    {'musculoskeletal_states'}, {'moment_arm'});





                       %% 2. De Groote equation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Groote F De, Kinney AL, Rao A V, Fregly BJ. Evaluation of Direct Collocation
% Optimal Control Problem Formulations for Solving the Muscle Redundancy Problem.
% Ann Biomed Eng. 2016; 44: 2922–2936. https://doi.org/10.1007/s10439-016-1591-9
% PMID: 27001399

%%                          single muscle 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 2.1 Parameters and input of the equation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 2.1.1 Muscle Tendon Parameters (ℓom, φo, Fom, ℓst) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
optimalFiberLength = SX.sym('Optimal_fiber_length', nMuscles);
phi0 = SX.sym('Pennation_angle_at_muscle_optimal_fiber_length', nMuscles);
maximalIsometricForce = SX.sym('Maximal_isometric_muscle_force', nMuscles);
tendonSlackLength = SX.sym('Tendon_slack_length', nMuscles);

muscleTendonParameters =  vertcat(optimalFiberLength, phi0, maximalIsometricForce, tendonSlackLength) ; % stock muscleTendonParameters

    % 2.1.2 Muscle tendon states input of the muscle tendon equation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a = SX.sym('Muscle_activation', nMuscles);
LUMT = SX.sym('UMT_length', nMuscles);

    % 2.1.3 rooted variables  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fiberLength = SX.sym('Fiber_length', nMuscles);
tendonLength = SX.sym('Tendon_length', nMuscles);
pennationAngle = SX.sym('Tendon_length', nMuscles);

rootedvariables = vertcat(fiberLength, tendonLength,pennationAngle);

    %% 2.2 Activation Dynamics (Not use in our modelisation) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% muscle activation is described by two nonlinear, first ODE
% parameter value
% taua = 0.015
% taud = 0.060
% b = 0.1
% taua = 0.015; %SX.sym('Activation_time_constant',nMuscles);
% taud = 0.060; %SX.sym('Desactivation_time_constant',nMuscles);
% b = 0.1; %SX.sym('Transition_smoothness', nMuscles);
% 
% fa = 0.5 .* tanh(b .* (e-a)) ;
% da_dt = ((1 ./ taua .* (0.5 + 1.5 .* a)) .* (fa + 0.5) + ...
%     ((0.5 + 1.5 .* a) ./ (taud)) .* (-fa + 0.5)) .* (e - a);
% ode = struct('x',a, 'u', e, 'ode', da_dt);
% tgrid = linspace(0, 1, 100);
% activation_dynamics = integrator('activation_dynamics', 'cvodes', ode, tgrid(1), tgrid(1:end));
% % example: 
% next_a = activation_dynamics('x0',[0;0;0],'u',[0.5; 0.7; 1]);
%
    %% 2.3 Muscle-Tendon Architecture Equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% muscle tendon architectural parameter
% tendon length --> input in function of tendon force 
% fiber length --> input in funcction of muscle force 
% muscle length --> input in muscle tendon equilibrium => length tnendon =
% length muscle = length fiber * cos(pennation angle) 

normalizedFiberLength = fiberLength ./ optimalFiberLength; 
normalizedTendonLength = tendonLength ./ tendonSlackLength; 

% muscleLength = umtLength - tendonLength; % fiberLength .* cos(pennationAngle); 

    %% 2.4  Muscle-Tendon Forces Equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2.4.2 Tendon Force Equation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tendon force-length (S1)
%%%%%%%%%%%%%%%%%%%%%%%%%%
kT = 35; c1 = 0.200; c2 = 0.995; c3 = 0.250; % tendon parameters
c3 = c3 - 0.0117508; % to avoid negative value 

normalizedTendonForcePart1 =  0 ;
normalizedTendonForcePart2 = (c1 .* exp(kT .* (normalizedTendonLength - c2)) - c3) ;
normalizedTendonForce = if_else(normalizedTendonLength < 1, ...
    normalizedTendonForcePart1, ...
    normalizedTendonForcePart2); % if normalized length under 0 the force = 0            % Normalized equation 

tendonForce = normalizedTendonForce.* maximalIsometricForce ; % Non-normalized equation

%% 2.4.2 Muscle Forces Equations 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Active force force-length (S2) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

b11 = 0.814483478343008 ; b21 = 1.055033428970575 ; b31 = 0.162384573599574 ; b41 = 0.063303448465465 ; % first Gaussian coefficents
b12 = 0.433004984392647 ; b22 = 0.716775413397760; b32 = -0.029947116970696 ; b42 = 0.200356847296188 ; % second Gaussian coefficents
b13 = 0.100 ; b23 = 1.000 ; b33 = 0.5 * sqrt(0.5) ; b43 = 0.000 ; % third Gaussian coefficents

num3 = normalizedFiberLength - b23;
den3 = b33 + b43 * normalizedFiberLength;
FMtilde3 = b13 * exp(-0.5 * num3.^2 ./ den3.^2);
num1 = normalizedFiberLength - b21;
den1 = b31 + b41 * normalizedFiberLength;
FMtilde1 = b11 * exp(-0.5 * num1.^2 ./ den1.^2);
num2 = normalizedFiberLength - b22;
den2 = b32 + b42 * normalizedFiberLength;
FMtilde2 = b12 * exp(-0.5 * num2.^2 ./ den2.^2);

normalizedMuscleActiveForceLength = FMtilde1 + FMtilde2 + FMtilde3; % Normalized equation 
muscleActiveForce = a .* normalizedMuscleActiveForceLength .* maximalIsometricForce ; % Non-normalized equation

% Passive muscle force-length (S3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kpe = 4.0 ; e0 = 0.6 ; % muscle parameters (oppti?) 

normalizedMusclePassiveForcePart1 =  0 ;
normalizedMusclePassiveForcePart2 = (exp(((kpe .* (normalizedFiberLength - 1))./e0)) - 1)./ (exp(kpe) - 1) ;
normalizedMusclePassiveForce = if_else(normalizedFiberLength < 1, ...
    normalizedMusclePassiveForcePart1, ...
    normalizedMusclePassiveForcePart2); % if normalized length under 0 the force = 0            % Normalized equation 

musclePassiveForce = normalizedMusclePassiveForce .* maximalIsometricForce ; % Non-normalized equation

% Passive muscle force-velocity (S4)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

normalizedMuscleForceVelocity = 1 ; % vitesse = 0

% Passive and active muscle force-length (S2, S3, S4)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

normalizedMuscleForce = a .* normalizedMuscleActiveForceLength .* normalizedMuscleForceVelocity + normalizedMusclePassiveForce ;   
muscleForce =  normalizedMuscleForce .* maximalIsometricForce ; 

% casadi fuction for vizualisation (single muscle)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (S1)
getNormalizedTendonForce = Function('getNormalizedTendonForce', ...
    {tendonLength(1), tendonSlackLength(1)}, {normalizedTendonForce(1)}, ...
    {'Tendon_length','Tendon_slack_length'}, {'normalizedTendonForce'}) ;
% (S2)
getNormalizedMuscleActiveForce = Function('getNormalizedMuscleActiveForce', ...
    {fiberLength(1), optimalFiberLength(1)}, {normalizedMuscleActiveForceLength(1)}, ...
    {'Fiber_length','Optimal_fiber_length'}, {'normalizedMuscleActiveForceLength'}) ;
% (S3)
getNormalizedMusclePassiveForce = Function('getNormalizedMusclePassiveForce', ...
    {fiberLength(1), optimalFiberLength(1)}, {normalizedMusclePassiveForce(1)}, ...
    {'Fiber_length','Optimal_fiber_length'}, {'normalizedMusclePassiveForce'}) ;
% (S2, S3, S4)
getNormalizedMuscleForce = Function('getNormalizedMuscleForce', ...
    {a(1),fiberLength(1), optimalFiberLength(1)}, {normalizedMuscleForce(1)}, ...
    {'Muscle_activation','Fiber_length','Optimal_fiber_length'}, {'normalizedMuscleForce'}) ;

                       %% 3. equilibrium functions 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %% 3.1  inputed and rooted variables of equilibrium functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% equilibrium function (for all muscles)
g5 = LUMT - ((cos(pennationAngle) .* fiberLength) + tendonLength) ;        % length equilibrium 
g6 = (optimalFiberLength .* sin(phi0)) - (fiberLength .* sin(pennationAngle)) ; % deforamtion 
g7 = muscleForce .* cos(pennationAngle) - tendonForce ;                    % forces equilibrium

unknown  = vertcat(tendonLength, fiberLength, pennationAngle) ;            % what we want to get 
known = vertcat(a, LUMT, muscleTendonParameters) ;                         % what we know 

% option of the solver 
opts_kinsol = struct("constraints", ones(9, 1), ...    
    'abstol', 1e-12, ...
    'u_scale', 1./[0.1, 0.01, 0.001, 0.1, 0.01, 0.001, 0.1, 0.01, 0.001]',...
    'max_iter',1000,...
    'error_on_fail', false,...
    'disable_internal_warnings', true);  %,...'iterative_solver', 'bcgstab', ...%'bcgstab', ...
    % 'print_level',3); % 1 means >= 0 
    
equilibriumError = Function('equilibriumError', {unknown, known}, {vertcat(g5, g6, g7)},{'x', 'p'}, {'residuals'}); 
equilibrateMuscleTendon = rootfinder('equilibrateMuscleTendon','kinsol',equilibriumError, opts_kinsol) ;

% single muscle 
% option of the solver 
opts_kinsol_single_muscle = struct("constraints", ones(3, 1), ...    
    'abstol', 1e-15, ...
    'u_scale', 1./[0.1, 0.04, 0.1]',...
    'max_iter',10000,...
    'error_on_fail', false,...
    'disable_internal_warnings', true); %,...  %,...'iterative_solver', 'bcgstab', ...%'bcgstab', ...
    %'print_level',3); % 1 means >= 0 

unknown  = vertcat(tendonLength(1), fiberLength(1), pennationAngle(1)) ;   % what we want to get 
known = vertcat(a(1), LUMT(1), muscleTendonParameters(1),muscleTendonParameters(4),muscleTendonParameters(7),muscleTendonParameters(10)) ;                % what we know 

equilibriumErrorSingleMuscle = Function('equilibriumErrorSingleMuscle', {unknown, known}, {vertcat(g5(1), g6(1), g7(1))},{'x', 'p'}, {'residuals'}); 
equilibrateMuscleTendonSingleMuscle = rootfinder('equilibrateMuscleTendonSingleMuscle','kinsol',equilibriumErrorSingleMuscle, opts_kinsol_single_muscle) ;


%% 2.5  Casadi functions about model Muscle-Tendon Forces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
neuromusculoskeletal_state = vertcat(a, q, skeleton) ; 

all_states = vertcat(neuromusculoskeletal_state,rootedvariables);

% about tendon forces
getTendonForce = Function('getTendonForce', ...
    {all_states, muscleTendonParameters}, {tendonForce}, ...
    {'all_states','muscle_tendon_parameters'}, {'tendon_force'}) ;

getNormalizeTendonForce = Function('getNormalizeTendonForce', ...
    {all_states, muscleTendonParameters}, {normalizedTendonForce}, ...
    {'all_states', 'muscleTendonParameters'}, {'NormalizedTendonForce'});

% about muscle forces 
getMuscleForce = Function('getMuscleForce', ...
    {all_states, muscleTendonParameters}, {muscleForce}, ...
    {'all_states','muscle_tendon_parameters'}, {'muscle_force'}) ;

getMusclePassiveForce = Function('getMusclePassiveForce', ...
    {all_states, muscleTendonParameters}, {musclePassiveForce}, ...
    {'all_states','muscle_tendon_parameters'}, {'MusclePassiveForce'}) ;

getMuscleActiveForce = Function('getMuscleActiveForce', ...
    {all_states, muscleTendonParameters}, {muscleActiveForce}, ...
    {'all_states','muscle_tendon_parameters'}, {'muscleActiveForce'}) ;

% about lenthening 
normalizeTendonLength = Function('normalizeTendonLength', ...
    {all_states, muscleTendonParameters}, {normalizedTendonLength}, ...
    {'all_states','muscle_tendon_parameters'}, {'length_tendon'});

normalizeFiberLength = Function('normalizeFiberLength', ...
    {all_states, muscleTendonParameters}, {normalizedFiberLength}, ...
    {'all_states','muscle_tendon_parameters'}, {'normalizeFiberLength'});



                       %% 4. Computing Joint Moments and Angles 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%  Moment articulaire Mj(θ, t) = ∑ i=1 m (ri(θ) ⋅ Fimt(θ, t)) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
jointMoment = momentArm' * tendonForce ;
jointMoment = jointMoment(end,end) ;
getJointMoment = Function('getJointMoment', ...
    {all_states, muscleTendonParameters}, {jointMoment}, ...
    {'all_states', 'muscle_tendon_parameters'}, {'jointMoment'});

FT = SX.sym('tendon_Forces', nMuscles);

jointMoment = momentArm' * FT ;
jointMoment = jointMoment(end,end) ;
getJointMoment2 = Function('getJointMoment', ...
    {musculoskeletal_states, FT}, {jointMoment}, ...
    {'musculoskeletal_states', 'FT'}, {'jointMoment'});







%% 5. struct that stock casadi function + Definition 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
casadiFun = struct(...
    'ForwardKinematics',ForwardKinematics,...                              % about musculoskeletal model (spatial configuration)
    'getUMTLength', getUMTLength,...
    'getMomentArm', getMomentArm,...                                       
    'getTendonForce', getTendonForce, ...                                  % about forces
    'getMuscleForce', getMuscleForce,...
    'getMusclePassiveForce',getMusclePassiveForce,...
    'getMuscleActiveForce', getMuscleActiveForce,...
    'getNormalizeTendonForce', getNormalizeTendonForce,...                
    'normalizeTendonLength',normalizeTendonLength,...                      % about lengthening
    'normalizeFiberLength', normalizeFiberLength, ...                      
    'equilibriumError',equilibriumError,...                                % about equilibrium
    'equilibrateMuscleTendon', equilibrateMuscleTendon,...
    'equilibriumErrorSingleMuscle',equilibriumErrorSingleMuscle,...        
    'equilibrateMuscleTendonSingleMuscle', equilibrateMuscleTendonSingleMuscle, ... 
    'getJointMoment',getJointMoment,...
    'getJointMoment2',getJointMoment2);                                    % about torque 

% function for vizualization
vizualizationFun = struct('getNormalizedTendonForce',getNormalizedTendonForce,...
    'getNormalizedMusclePassiveForce', getNormalizedMusclePassiveForce,...
    'getNormalizedMuscleActiveForce', getNormalizedMuscleActiveForce,...                                       
    'getNormalizedMuscleForce', getNormalizedMuscleForce); 

definition = ["a : Neuromuscular activation (between 0 and 1)";...
    "q : Spatial skeleton configuration (x, y, z, theta_hip, theta_knee, theta_ankle)" ;...
    "skeleton : musculoskeletal parameters (segment_geometry, muscle_insertion, Local_ViaPoint_tibialis_anterior)";...
    "neuromusculoskeletal_state : neuromusculoskeletal parameters(a, q, known_parameters)" ;...
    "muscleTendonParameters : Muscle Tendon Parameters (ℓom, φo, Fom, ℓst) " ;...
    "p : neuromusculoskeletal parameters + Muscle Tendon Parameters (neuromusculoskeletal_state, muscleTendonParameters) " ] ; 

%% Unknown parameters for NLP
unknown_parameters = horzcat(optimalFiberLength, phi0, maximalIsometricForce, tendonSlackLength);
end

function [rotation_matrix] = Rototranslation_Rz(x,y,z, rotZ)
%import casadi.SX
import casadi.*

rotation_matrix = SX.eye(4);
rotation_matrix(1:2, 1:2) = [ cos(rotZ), -sin(rotZ) ; ...
    sin(rotZ), cos(rotZ)] ;
rotation_matrix(1:3,4) = [x, y, z]';
end

function [point_in0] = Rototranslate(R1to0, point_in1)
point_in0 = R1to0 * [point_in1;1];
point_in0 = point_in0(1:3);
end
