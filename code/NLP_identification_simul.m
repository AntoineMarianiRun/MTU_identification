% import casadi
[ret, name] = system('hostname');
if strcmp(name, '942-27984')
    addpath('C:\Users\Stage\Desktop\Neuromusculoskeletal Modeling\Casadi')
elseif strcmp(name, 'xxx')
    addpath('/Users/mickaelbegon/Downloads/casadi-3.6.3-osx64-matlab2018b/')
end
addpath('C:\Users\Stage\Desktop\Neuromusculoskeletal Modeling\Casadi')

import casadi.*


% Call data_generator
hypothetical_data_generator ;

% Add noise
%%

% Start with an empty NLP
w={}; %variables
w0 = []; % initial guess
lbw = []; % lower bound of variable
ubw = []; % upper bound of variable
J = 0; % cost | objective function
gg={}; % constraints
lbg = []; % lower bound of constraints
ubg = []; % upper bound of constraints
gg0 = [] ; 
UP = horzcat(unknown_parameters(:,1)',unknown_parameters(:,2)',...
    unknown_parameters(:,3)',unknown_parameters(:,4)') ;

% create muscle parameters variables
% w = { w{:}, UP};
X = {} ; 

% Weightings in cost function
W_torque = 1; % 1 Nm
W_length = 0.005; %todo 1Nm correspond à 50 mm
W_angle = (3/180) * pi;  % 1 Nm correspond à 3 deg


for trial = 1 : ntrials % for 1 to nb trials

    % create variable for each trial
    w = { w{:}, UP};
    w0 =  [w0, muscle_tendon_parameters_num ]; 
    lbw = [lbw, muscle_tendon_parameters_num * 0.5]; % lower bound of variable
    ubw = [ubw, muscle_tendon_parameters_num * 2]; % upper bound of variable
    lbg = [lbg, zeros(1,6)]; % lower bound of constraints
    ubg = [ubg, [0.05, 0.05, 0.05 ,muscle_tendon_parameters_num(1:3)*1.8]]; % upper bound of constraints
    gg0 = [gg0, 0.01, 0.01, 0.01 ,muscle_tendon_parameters_num(1:3)] ; % initial guess

    fiberLength_k = SX.sym(['Fiber_length_' num2str(trial)], nMuscles);
    tendonLengthening_k = SX.sym(['Tendon_Lengthening_' num2str(trial)], nMuscles);
    x = vertcat((tendonLengthening_k),(fiberLength_k)) ; 

    X = {X{:},x'};


   
    data = Data(trial,:);

        % known variables 
    a_trial = data(4:6) ; % muscle activation during the trial
    q_trial = [0, 0, 0, 0, data(2:3)] ; % skeleton configuration during the trial
    musculoskeletal_states_trial = [q_trial, known_parameters_num ] ; % muscleskeleton configuration during the trial
    neuromusculoskeletal_states_trial = [a_trial, musculoskeletal_states_trial] ; % Neuromusculoskeletal states
    p_trial = horzcat(neuromusculoskeletal_states_trial, w{trial}) ; % states of 'rooted equation 'g'
    umtLength = casadiFun.getUMTLength(musculoskeletal_states_trial) ; % umt length 

    % constraints
    [g0,g1]  =  muscletendonequation(a_trial,fiberLength_k,tendonLengthening_k,unknown_parameters,umtLength) ; 
    g = Function('g', {x, p}, {vertcat(g0,g1)},{'x', 'p'}, {'residuals'}) ;
    contraints = g(x,p_trial);

    % equation of estimated variables (Torque, fiber length, pennation angle)
    tendonLengthening_trial = contraints(1:nMuscles)' ; 
    fiberLength_trial = contraints(nMuscles+1:end)' ;
    rootedvariables_trial = [fiberLength_trial ,tendonLengthening_trial] ;
    all_states_trial = [neuromusculoskeletal_states_trial, rootedvariables_trial] ;

    temp = casadiFun.getJointMoment(all_states_trial,  w{trial}) ;
    Torque_simulated = temp(end) ; % Torque equation 
    FiberLength_simulated = fiberLength_trial ; % fiber length equation
    phi_simulated = casadiFun.getPennationAngle(all_states_trial,w{trial})' ; % pennation angle equation

    % objective
     J = J + W_torque * (data(1) - Torque_simulated)^2; %add error on joint torque
     J = J + W_length * sum((data(7:9) - FiberLength_simulated).^2);% add error on tendon length 
     J = J + W_angle * sum((data(10:12) - phi_simulated).^2); %add error pennation angle

     gg = { gg{:}, contraints}; % contraints functions 
end

%%
prob = struct('x', [horzcat(w{1}),horzcat(X{:})], 'f', J , 'g',vertcat(gg{:})); 
solver = nlpsol('solver', 'ipopt', prob);

 W0 = [w0(1:12), gg0] ; % inital guess vector 
sol = solver('x0', W0, 'lbx', [lbw(1,1:12), lbg], 'ubx', [ubw(1,1:12),ubg] , ...
    'lbg', lbg, 'ubg', ubg);

w_opt = full(sol.x);
nparam = 12 ; 
param_opt = w_opt(1:nparam);
difference = w0(1:12) - param_opt ; 
