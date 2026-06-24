import os
import import_functions
import manipfun
import useful
import matplotlib.pyplot as plt
import numpy as np
import model_plot as pm_

#############################################################################
#     1. Organized path and files in a dictionary
#############################################################################
current_folder =  os.getcwd()
measured_data_folder = os.path.join(current_folder, "num_data")
subject_name = 'BOC_09'
condition_name = 'sDyn' # 'sDyn' or 'sPra'
trial_name = subject_name + '_' + condition_name

folders = {
    "main":current_folder,
    "fun": os.path.join(current_folder, "pyFun"),
    "osim_model": os.path.join(current_folder, "osimModel"),
    "measured_data": os.path.join(measured_data_folder, subject_name),
    "sim_data": os.path.join(current_folder, "simulated_data"),
}

#############################################################################
#     2. neuro-musculo-skeletal: generic model
#############################################################################
"""
        # 2.1 neuro-musculo-skeletal (ℓom, φo, Fom, ℓst)
# 2.1.1 Import from Opensim Musculoskeletal Geometry and Generic Muscle Tendon Parameters
mtu_params=['l0m', 'phi0', 'f0m', 'lst']


# 2.1.2 Scrape from .osim file geometry of bodies and muscle-tendon parameters
skeleton_num, muscle_tendon_parameters_num = import_functions.get_model_osim_generic(folders["osim_model"], "wholebody.osim",mtu_params = mtu_params)


# 2.1.3 Get kinematic and dynamic equations (muscle tendon equation from De
# Groote) --> Hill type model --> Fmt = f(a, ℓmt, νmt; Fom, ℓom, ℓst, φo).
# Note that in our model we ignore :
#       - fiber contraction velocity (νmt = 1)
#       - and electromechanical delay (a(t) = e(t))
casadi_function, unknown_parameters, definition = useful.get_model_equation()


# 2.1.4 Test the model
# 2.1.4.1 skeletal state (q)
q1 = 0 # x coord pelvis
q2 = 0 # y coord pelvis
q3 = 0 # z coord pelvis
q4 = 0 # phi angle pelvis
q5 = 0 # phi angle knee
q6 = -34.3 # phi angle ankle

q4 = (q4 / 180) * np.pi
q5 = (q5 / 180) * np.pi
q6 = (q6 / 180) * np.pi

# 2.1.4.2 concatenate
q_num = [q1,q2,q3,q4,q5,q6]

# 2.1.4.3 muscle state (a)
a_num = [0.0, # activation tibialis
         0.0, # activation soleus
         0.0] # activation gast

# 2.1.4.4 test
useful.test_model(skeleton_num,muscle_tendon_parameters_num,casadi_function,a_num,q_num)


# 2.1.4 plot and verification of the muscle and tendon dynamic equation
fig, axes = useful.plot_force_length_activation_3d(
    a=[0.7, 0.7, 0.5],
    fiber_length=[0.098*1.3, 0.08, 0.12],
    tendon_length=[0.223*1.04, 0.18, 0.25],
    muscle_tendon_parameters=muscle_tendon_parameters_num,           # vecteur de 12
    get_fiber_force_from_fiber_length=casadi_function["get_fiber_force_from_fiber_length"],
    get_tendon_force_from_tendon_length=casadi_function["get_tendon_force_from_tendon_length"],
    muscle_idx=1,
    muscle_names=["tibialis", "Soleus", "Gastroc"],
    normalize=True,
)

plt.savefig("force_length_activation_3d.pdf", dpi=300, bbox_inches="tight")

# 2.1.5 interactive model test
useful.interactive_model(skeleton_num, muscle_tendon_parameters_num, casadi_function)


"""
        # 2.2 neuro-musculo-skeletal (ℓom, φo, Fom, kpe, ℓst, kt)
# 2.2.1 Import from Opensim Musculoskeletal Geometry and Generic Muscle Tendon Parameters
mtu_params=['l0m', 'phi0', 'f0m', 'km', 'lst','kt']


# 2.2.2 Scrape from .osim file geometry of bodies and muscle-tendon parameters
skeleton_num, muscle_tendon_parameters_num = import_functions.get_model_osim_generic(folders["osim_model"], "wholebody.osim",mtu_params = mtu_params)


# 2.2.3 Get kinematic and dynamic equations (muscle tendon equation from De
# Groote) --> Hill type model --> Fmt = f(a, ℓmt, νmt; [ℓom, φo, Fom, kpe, ℓst, kt]).
# Note that in our model we ignore :
#       - fiber contraction velocity (νmt = 1)
#       - and electromechanical delay (a(t) = e(t))
param_config = {
    'l0m': 'sym',
    'phi0': 'sym',
    'f0m': 'sym',
    'km': 'sym',
    'lst': 'sym',
    'kt': 'sym',
}
casadi_function, unknown_parameters, definition = useful.get_model_equation(param_config)


# 2.2.4 Test the model
# 2.2.4.1 skeletal state (q)
q1 = 0 # x coord pelvis
q2 = 0 # y coord pelvis
q3 = 0 # z coord pelvis
q4 = 0 # phi angle pelvis
q5 = 80 # phi angle knee
q6 = 30 # phi angle ankle
q4 = (q4 / 180) * np.pi
q5 = (q5 / 180) * np.pi
q6 = (q6 / 180) * np.pi

# 2.2.4.2 concatenate
q_num = [q1,q2,q3,q4,q5,q6]

# 2.2.4.3 muscle state (a)
a_num = [0.2, # activation tibialis
         0.2, # activation soleus
         0.2] # activation gast


# ----------------------------------------------------------------------
# Exemple d'utilisation
# ----------------------------------------------------------------------
    # ordre = [p for p in mtu_params if sym] = l0m, phi0, f0m, km, lst, kt

# ----------------------------------------------------------------------
# Exemple (stubs numpy imitant la signature CasADi : entrées[3], params[18])
# ----------------------------------------------------------------------

a0 = np.array([0.3, 0.5, 0.2])
fl0 = np.array([0.10, 0.13, 0.085])
tl0 = np.array([0.205, 0.182, 0.255])


w = pm_.force_length_sliders(
    param_config,
    muscle_tendon_parameters_num,
    a0, fl0, tl0,
    casadi_function['get_fiber_force_from_fiber_length'],
    casadi_function['get_tendon_force_from_tendon_length'],
    muscle_idx=1,
    muscle_names=["Tibialis", "Soleus", "Gastrocnemius"],
    param_layout="by_type",
)


"""
# 2.2.4.4 test
pm_.interactive_model(skeleton_num, muscle_tendon_parameters_num, casadi_function)


# 2.2.5 interactive model test
useful.interactive_model(skeleton_num, muscle_tendon_parameters_num, casadi_function)

"""
#############################################################################
#     3. neuro-musculo-skeletal: scaled model
#############################################################################
"""
        # 3.1 neuro-musculo-skeletal (ℓom, φo, Fom, ℓst)
# 3.1.1 Import from Opensim Musculoskeletal Geometry and Generic Muscle Tendon Parameters
mtu_params=['l0m', 'phi0', 'f0m', 'lst']


# 3.1.2 Scrape from .osim file geometry of bodies and muscle-tendon parameters
skeleton_num, muscle_tendon_parameters_num = import_functions.get_model_osim_scaled(folders["measured_data"], f"{subject_name}.osim",mtu_params = mtu_params)


# 3.1.3 Get kinematic and dynamic equations (muscle tendon equation from De
# Groote) --> Hill type model --> Fmt = f(a, ℓmt, νmt; Fom, ℓom, ℓst, φo).
# Note that in our model we ignore :
#       - fiber contraction velocity (νmt = 1)
#       - and electromechanical delay (a(t) = e(t))
casadi_function, unknown_parameters, definition = useful.get_model_equation()


# 3.1.4 Test the model
# 3.1.4.1 skeletal state (q)
q1 = 0 # x coord pelvis
q2 = 0 # y coord pelvis
q3 = 0 # z coord pelvis
q4 = 0 # phi angle pelvis
q5 = 80 # phi angle knee
q6 = 30 # phi angle ankle

q4 = (q4 / 180) * np.pi
q5 = (q5 / 180) * np.pi
q6 = (q6 / 180) * np.pi

# 3.1.4.2 concatenate
q_num = [q1,q2,q3,q4,q5,q6]

# 3.1.4.3 muscle state (a)
a_num = [0.2, # activation tibialis
         0.2, # activation soleus
         0.2] # activation gast

# 3.1.4.4 test
useful.test_model(skeleton_num,muscle_tendon_parameters_num,casadi_function,a_num,q_num)


# 3.1.5 interactive model test
useful.interactive_model(skeleton_num, muscle_tendon_parameters_num, casadi_function)


        # 3.2 neuro-musculo-skeletal (ℓom, φo, Fom, kpe, ℓst, kt)
# 3.2.1 Import from Opensim Musculoskeletal Geometry and Generic Muscle Tendon Parameters
mtu_params=['l0m', 'phi0', 'f0m', 'km', 'lst','kt']


# 3.2.2 Scrape from .osim file geometry of bodies and muscle-tendon parameters
skeleton_num, muscle_tendon_parameters_num = import_functions.get_model_osim_scaled(folders["measured_data"], f"{subject_name}.osim",mtu_params = mtu_params)


# 3.2.3 Get kinematic and dynamic equations (muscle tendon equation from De
# Groote) --> Hill type model --> Fmt = f(a, ℓmt, νmt; [ℓom, φo, Fom, kpe, ℓst, kt]).
# Note that in our model we ignore :
#       - fiber contraction velocity (νmt = 1)
#       - and electromechanical delay (a(t) = e(t))
param_config = {
    'l0m': 'sym',
    'phi0': 'sym',
    'f0m': 'sym',
    'km': 'sym',
    'lst': 'sym',
    'kt': 'sym',
}
casadi_function, unknown_parameters, definition = useful.get_model_equation(param_config)


# 3.2.4 Test the model
# 3.2.4.1 skeletal state (q)
q1 = 0 # x coord pelvis
q2 = 0 # y coord pelvis
q3 = 0 # z coord pelvis
q4 = 0 # phi angle pelvis
q5 = 80 # phi angle knee
q6 = 30 # phi angle ankle
q4 = (q4 / 180) * np.pi
q5 = (q5 / 180) * np.pi
q6 = (q6 / 180) * np.pi

# 3.2.4.2 concatenate
q_num = [q1,q2,q3,q4,q5,q6]

# 3.2.4.3 muscle state (a)
a_num = [0.2, # activation tibialis
         0.2, # activation soleus
         0.2] # activation gast

# 3.2.4.4 test
useful.test_model(skeleton_num,muscle_tendon_parameters_num,casadi_function,a_num,q_num)


# 3.2.5 interactive model test
useful.interactive_model(skeleton_num, muscle_tendon_parameters_num, casadi_function)

"""
#############################################################################
#    4. hypothetical data generator
#############################################################################
"""
# 4.1.1 Import from Opensim Musculoskeletal Geometry and Generic Muscle Tendon Parameters
mtu_params=['l0m', 'phi0', 'f0m', 'lst']


# 4.1.2 Scrape from .osim file geometry of bodies and muscle-tendon parameters
skeleton_num, muscle_tendon_parameters_num = import_functions.get_model_osim_scaled(folders["measured_data"], f"{subject_name}.osim",mtu_params = mtu_params)


# 4.1.3 Get kinematic and dynamic equations (muscle tendon equation from De
# Groote) --> Hill type model --> Fmt = f(a, ℓmt, νmt; Fom, ℓom, ℓst, φo).
# Note that in our model we ignore :
#       - fiber contraction velocity (νmt = 1)
#       - and electromechanical delay (a(t) = e(t))
casadi_function, unknown_parameters, definition = useful.get_model_equation()


# 4.2 data generation
header, hypothetical_data = useful.hypothetical_data_generator(skeleton_num, muscle_tendon_parameters_num, casadi_function)


# 4.3 visualization
manipfun.plot_data(hypothetical_data, muscle_names=['tibialis', 'soleus', 'gastrocnemius'])


# 4.4 save data in an Excel file
useful.save_data_to_xlsx(hypothetical_data,header,folders["sim_data"],'dataset_hypothetical')


# 4.5 load data from an Excel file
data = useful.load_data_from_xlsx(folders["sim_data"], 'dataset_hypothetical', header)


# 4.6 add tendon length to data set and import test data set
data = manipfun.add_tendon_length_to_data(data,skeleton_num, casadi_function)


# 4.6 comparison
stats = useful.compare_datasets(
    hypothetical_data, data, header,
    label_a='original', label_b='reloaded'
)
"""
#############################################################################
#    5. NLP  NonLinear Programming optimisation problem (ℓom, φo, Fom, ℓst)
#                           NUMERIC VALIDATION
#############################################################################
"""
# 5.1. import the model
# 5.1.1 Import from Opensim Musculoskeletal Geometry and Generic Muscle Tendon Parameters
mtu_params=['l0m', 'phi0', 'f0m', 'lst','kt']

param_config = {
    'l0m': 'sym',
    'phi0': 'sym',
    'f0m': 'sym',
    'km': 'fixed',
    'lst': 'sym',
    'kt': 'sym',
}

# 5.1.2 Scrape from .osim file geometry of bodies and muscle-tendon parameters
skeleton_num, muscle_tendon_parameters_num = import_functions.get_model_osim_scaled(folders["measured_data"], f"{subject_name}.osim",mtu_params = mtu_params)


# 5.1.3 Get kinematic and dynamic equations (muscle tendon equation from De
# Groote) --> Hill type model --> Fmt = f(a, ℓmt, νmt; Fom, ℓom, ℓst, φo).
# Note that in our model we ignore :
#       - fiber contraction velocity (νmt = 1)
#       - and electromechanical delay (a(t) = e(t))
casadi_function, unknown_parameters, definition = useful.get_model_equation(param_config=param_config)


# 5.2 get data
# 5.2.1 load data from an Excel file
header = [
    'ankle_torque',
    'q_knee', 'q_ankle',
    'a_tibialis', 'a_soleus', 'a_gastrocnemius',
    'fiber_length_tibialis', 'fiber_length_soleus', 'fiber_length_gastrocnemius',
    'pennation_angle_tibialis', 'pennation_angle_soleus', 'pennation_angle_gastrocnemius',
    'tendon_length_tibialis', 'tendon_length_soleus', 'tendon_length_gastrocnemius'
]

hypothetical_data = useful.load_data_from_xlsx(folders["sim_data"], 'dataset_hypothetical', header)


# 5.2.2 add tendon length to data set and import test data set
hypothetical_data = manipfun.add_tendon_length_to_data(hypothetical_data,skeleton_num, casadi_function)

# 5.2.3 plot datas (verif)
manipfun.plot_data(hypothetical_data, muscle_names=['tibialis', 'soleus', 'gastrocnemius'])


#  5.3 set the initial guess perfect
initial_guess, upper_band, lower_band, param_index = useful.get_initial_guess(
    muscle_tendon_parameters_num, hypothetical_data, param_config,
    'scaled', verbose=True)



#  5.4 optimisation problem with perfect initial guess
muscle_tendon_parameters_opt_xO_parf = useful.optimization_nlp(
    hypothetical_data,
    muscle_tendon_parameters_num,
    lower_band,
    upper_band,
    skeleton_num,
    muscle_tendon_parameters_num,
    unknown_parameters,
    casadi_function,
    param_index)


#  5.5 generate data with the optimized parameter
header, hypothetical_data_opt, path_csv = useful.generate_estimated_data(hypothetical_data, skeleton_num, muscle_tendon_parameters_num,
                            casadi_function,
                            output_dir='simulated_data', filename='data_estime_generic_with_perfect_guess.csv',
                            save_npy=False, save_csv=True, verbose=True)


stats = useful.compare_datasets(
    hypothetical_data, hypothetical_data_opt, header,
    label_a='generic', label_b='opti x0 perfect'
)


#  5.6 optimisation problem with initial guess according to our methods
initial_guess, upper_band, lower_band, param_index = useful.get_initial_guess(
    muscle_tendon_parameters_num, hypothetical_data, param_config,
    'measured', verbose=True)


muscle_tendon_parameters_opt_x0_rand = useful.optimization_nlp(
    hypothetical_data,
    initial_guess,
    lower_band,
    upper_band,
    skeleton_num,
    muscle_tendon_parameters_num,
    unknown_parameters,
    casadi_function,
    param_index)


muscle_tendon_parameters_opt_x0_rand = useful.optimization_nlp(
    hypothetical_data, initial_guess, lower_band, upper_band, skeleton_num,
    muscle_tendon_parameters_num, unknown_parameters, casadi_function,
    param_index)


#  5.7 generate data with the optimized parameter
header, hypothetical_data_opt, path_csv = useful.generate_estimated_data(hypothetical_data, skeleton_num, muscle_tendon_parameters_opt_x0_rand,
                            casadi_function,
                            output_dir='simulated_data', filename='data_estime_measured.csv',
                            save_npy=False, save_csv=True, verbose=True)

manipfun.plot_data(hypothetical_data_opt, muscle_names=['tibialis', 'soleus', 'gastrocnemius'])

stats = useful.compare_datasets(
    hypothetical_data, hypothetical_data_opt, header,
    label_a='original', label_b='opti'
)

rng = np.random.default_rng(seed=42)

hypothetical_data_noise = useful.add_noise(
    hypothetical_data,
    sigma_torque = 0.5,
    sigma_length=0.002,
    sigma_angle_deg=0.5,
    add_internal_noise=True,
    rng=rng,
)

hypothetical_data_noise = manipfun.add_tendon_length_to_data(hypothetical_data_noise,skeleton_num, casadi_function)

manipfun.plot_data(hypothetical_data_noise, muscle_names=['tibialis', 'soleus', 'gastrocnemius'])


muscle_tendon_parameters_opt_bruit = useful.optimization_nlp(
    hypothetical_data_noise, initial_guess, lower_band, upper_band, skeleton_num,
    muscle_tendon_parameters_num, unknown_parameters, casadi_function,
    param_index)

    
initial_guess, upper_band, lower_band = useful.get_initial_guess(muscle_tendon_parameters_num, hypothetical_data,
                                  strategy='measured')





# 5.8 conditioning and identifiability
initial_guess, upper_band, lower_band, param_index = useful.get_initial_guess(
    muscle_tendon_parameters_num, hypothetical_data, param_config,
    'scaled', verbose=True)

initial_guess = muscle_tendon_parameters_num

from monte_carlo_identification import main as run_mc

results_mc = run_mc(
    data_clean=hypothetical_data,
    initial_guess=initial_guess,
    lower_band=lower_band,
    upper_band=upper_band,
    skeleton_num=skeleton_num,
    muscle_tendon_parameters_num=muscle_tendon_parameters_num,
    unknown_parameters=unknown_parameters,
    casadi_function=casadi_function,
    optimization_nlp=useful.optimization_nlp,
    param_index=param_index,
    cfg={"n_mc": 100},
)


from run_identifiability import run_all

setup = dict(
    data=hypothetical_data,
    initial_guess=initial_guess,
    lower_band=lower_band,
    upper_band=upper_band,
    skeleton_num=skeleton_num,
    true_params=muscle_tendon_parameters_num,
    unknown_parameters=unknown_parameters,
    casadi_function=casadi_function,
    optimization_nlp=useful.optimization_nlp,
    param_index=param_index,
    param_config=param_config,
)

results_id = run_all(
    setup,
    analyses = ("profile"),
    profile_kwargs={"n_grid": 22, "rel_range": 0.2},
    sobol_kwargs={"n_base": 64},        # petit pour un premier essai rapide
)
 # analyses = ("conditioning", "profile", "sobol"),
"""
#############################################################################
#    6. NLP  NonLinear Programming optimisation problem (ℓom, φo, Fom, ℓst)
#                           EXPERIMENTAL VALIDATION
#                               train and test
#############################################################################

# 6.1 Import datasets
# 6.1.1 name of the measured data
header = [
    'ankle_torque',
    'q_knee', 'q_ankle',
    'a_tibialis', 'a_soleus', 'a_gastrocnemius',
    'fiber_length_tibialis', 'fiber_length_soleus', 'fiber_length_gastrocnemius',
    'pennation_angle_tibialis', 'pennation_angle_soleus', 'pennation_angle_gastrocnemius',
    'tendon_length_tibialis', 'tendon_length_soleus', 'tendon_length_gastrocnemius'
]

# 6.1.2 get data folder
osim_path, train_path, test_path = manipfun.select_folder_and_get_files() #  get the interest folder

# 6.1.3 name and path of the files
osim_folder, osim_name = manipfun.split_path_name(osim_path)
train_folder, train_name = manipfun.split_path_name(train_path)
test_folder, test_name = manipfun.split_path_name(test_path)

# 6.1.4 import the osim scaled as a python variable
mtu_params=['l0m', 'phi0', 'f0m', 'lst','kt']
param_config = {
    'l0m': 'sym',
    'phi0': 'sym',
    'f0m': 'sym',
    'km': 'fixed',
    'lst': 'sym',
    'kt': 'sym',
}

skeleton_num, muscle_tendon_parameters_num = import_functions.get_model_osim_scaled(osim_folder,osim_name,mtu_params = mtu_params)
casadi_function, unknown_parameters, definition = useful.get_model_equation(param_config=param_config)

# 6.1.5 training data set and import test data set
data_train = useful.load_data_from_xlsx(train_folder, train_name, header)
data_test = useful.load_data_from_xlsx(test_folder, test_name, header)

# 6.1.6 add tendon length to data set and import test data set
data_train = manipfun.add_tendon_length_to_data(data_train,skeleton_num, casadi_function)
data_test = manipfun.add_tendon_length_to_data(data_test,skeleton_num, casadi_function)

# 6.1.7 data visual verification

manipfun.plot_data(data_train, muscle_names=['tibialis', 'soleus', 'gastrocnemius'])
manipfun.plot_data(data_test, muscle_names=['tibialis', 'soleus', 'gastrocnemius'])



# 6.2 optimization
# 6.2.1 set the initial guess according to the test data

initial_guess, upper_band, lower_band, param_index = useful.get_initial_guess(
    muscle_tendon_parameters_num, data_train, param_config,
    'measured', verbose=True)


initial_guess, upper_band, lower_band, param_index = useful.get_initial_guess(
    muscle_tendon_parameters_num, data_train, param_config,
    'scaled', verbose=True)


# 4. Lancer la figure interactive

pm_.interactive_calibration_figure(
    data=data_test,
    skeleton_num=skeleton_num,
    casadi_function=casadi_function,
    param_index=param_index,
    param_initial=initial_guess,
    param_lower=lower_band,
    param_upper=upper_band,
    solve_all_muscles=useful.solve_all_muscles,
)

lower_band[12:15] *=0.01
initial_guess[12:15] *=1.3
initial_guess[9:11] *=1.3
upper_band[9:11] *= 2


# 6.2.2 optimisation problem

muscle_tendon_parameters_opt = useful.optimization_nlp(
    data_train,
    initial_guess,
    lower_band,
    upper_band,
    skeleton_num,
    muscle_tendon_parameters_num,
    unknown_parameters,
    casadi_function,
    param_index)


"""
muscle_tendon_parameters_opt = useful.optimization_nlp(
    data_train,
    initial_guess,
    lower_band,
    upper_band,
    skeleton_num,
    muscle_tendon_parameters_num,
    unknown_parameters,
    casadi_function,
    param_index)
    
import nlp_test

muscle_tendon_parameters_opt = nlp_test.optimization_nlp_raw_plain(
    data_train,
    initial_guess,
    lower_band,
    upper_band,
    skeleton_num,
    muscle_tendon_parameters_num,
    unknown_parameters,
    casadi_function,
    param_index)
"""
"""
useful.interactive_model(skeleton_num, muscle_tendon_parameters_opt, casadi_function)




# 6.2.3 opt parameters save
manipfun.save_muscle_tendon_parameters(muscle_tendon_parameters_opt,
                                  osim_folder,
                                  filename='muscle_tendon_parameters_saved',
                                  muscle_names=('tibialis', 'soleus', 'gastrocnemius'),
                                  save_csv=True,
                                  save_npy=True,
                                  verbose=True)

# 6.3 test
# 6.3.1 generate test data with the optimized parameter
header, data_opt, path_csv = useful.generate_estimated_data(data_test, skeleton_num, muscle_tendon_parameters_opt,
                            casadi_function,
                            output_dir='simulated_data', filename='data_estime_opt.csv',
                            save_npy=False, save_csv=False, verbose=True)



# 6.3.1 generate test data with the optimized parameter
header, data_opt, path_csv = useful.generate_estimated_data(data_test, skeleton_num, muscle_tendon_parameters_num,
                            casadi_function,
                            output_dir='simulated_data', filename='data_estime_opt.csv',
                            save_npy=False, save_csv=False, verbose=True)

# 6.3.2 visual verification of generate test data and test data
manipfun.plot_data(data_test, muscle_names=['tibialis', 'soleus', 'gastrocnemius'])
manipfun.plot_data(data_opt, muscle_names=['tibialis', 'soleus', 'gastrocnemius'])

# 6.3.2 differences between generate test data and measured test data
stats = useful.compare_datasets(
    data_test, data_opt, header,
    label_a='original', label_b='reloaded'
)


header, data_nopt, path_csv = useful.generate_estimated_data(data_test, skeleton_num, muscle_tendon_parameters_num,
                            casadi_function,
                            output_dir='simulated_data', filename='data_estime_generic.csv',
                            save_npy=False, save_csv=False, verbose=True)

manipfun.plot_data(data_test, muscle_names=['tibialis', 'soleus', 'gastrocnemius'])
manipfun.plot_data(data_nopt, muscle_names=['tibialis', 'soleus', 'gastrocnemius'])

stats = useful.compare_datasets(
    data_test, data_est_nopt, header,
    label_a='original', label_b='reloaded'
)



# 6.5 test the difference between measured data and data simulated with a generic model
# 6.5.1 data simulated with a generic model
header, simulated_data, path_csv = useful.generate_estimated_data(data_test, skeleton_num, muscle_tendon_parameters_num,
                            casadi_function,
                            output_dir='simulated_data', filename='data_estime_opt.csv',
                            save_npy=False, save_csv=False, verbose=True)

# 6.5.2 plot data simulated
manipfun.plot_data(simulated_data, muscle_names=['tibialis', 'soleus', 'gastrocnemius'])

# 6.5.3 stats
stats = useful.compare_datasets(
    data_test, simulated_data, header,
    label_a='original', label_b='reloaded'
)


mtu_params=['l0m', 'phi0', 'f0m','km', 'lst','kt']
skeleton_num, muscle_tendon_parameters_num = import_functions.get_model_osim_scaled(osim_folder,osim_name,mtu_params = mtu_params)
param_config = {
    'l0m': 'sym',
    'phi0': 'sym',
    'f0m': 'sym',
    'km': 'sym',
    'lst': 'sym',
    'kt': 'sym',
}
casadi_function, unknown_parameters, definition = useful.get_model_equation(param_config)


initial_guess = muscle_tendon_parameters_num
upper_band = muscle_tendon_parameters_num + muscle_tendon_parameters_num *.9
lower_band = muscle_tendon_parameters_num - muscle_tendon_parameters_num *.9


muscle_tendon_parameters_opt = useful.optimization_nlp(data_train,initial_guess,lower_band,upper_band,skeleton_num,muscle_tendon_parameters_num,unknown_parameters,casadi_function)

"""
#############################################################################
#    7. Simulation based on (skeleton_num, muscle_tendon_parameters_num, casadi_function,time_num, q5_num, q6_num, aTibialis_num, aSoleus_num, aGastrocnemius_num)
#############################################################################
"""
# find measured data
##########################
measured_data_folder_name = os.path.join(folders["measured_data"], trial_name + '.xlsx')
df = pd.read_excel(measured_data_folder_name) # import your file

# Assign columns to variables
time_series = df["time (s)"]      # pandas Series time
#q5_series = df["time (s)"]      # pandas Series knee
q6_series = df["ankle angle (deg)"]      # pandas Series ankle
a_tibialis_series = df["tibialis emg"]      # pandas Series ankle
a_soleus_series = df["soleus emg"]      # pandas Series ankle
a_gastrocnemius_series = df["gastrocnemius emg"]      # pandas Series ankle

# Example: Convert to Python lists
time_num = time_series.tolist()
q6_num = q6_series.tolist()
q6_num = [-x for x in q6_num] # change orientation frame
q5_num = [0] * len(q6_num)
a_tibialis_num = a_tibialis_series.tolist()
a_soleus_num = a_soleus_series.tolist()
a_gastrocnemius_num = a_gastrocnemius_series.tolist()


# change normalization of emg from % to ratio
# if any(x < 1 for x in a_tibialis_num):
    # divide every element by 100
a_tibialis_num = [x/100 for x in a_tibialis_num]

# if any(x < 1 for x in a_soleus_num):
    # divide every element by 100
a_soleus_num = [x/100 for x in a_soleus_num]

# if any(x < 1 for x in a_gastrocnemius_num):
    # divide every element by 100
a_gastrocnemius_num = [x / 100 for x in a_gastrocnemius_num]

# define activation limit

# simulation data
##########################
header_neuro_musculo_state, simulated_data_neuro_musculo_state, header_mtu_architecture, simulated_data_mtu_architecture, header_mtu_forces, simulated_data_mtu_forces = useful.simulation_(skeleton_num, muscle_tendon_parameters_num, casadi_function,time_num, q5_num, q6_num, a_tibialis_num, a_soleus_num, a_gastrocnemius_num)


# save simulated data
##########################
    # 1. Save neuro_musculo_state
    ##########################
# Save path
save_dir = os.path.join(folders["measured_data"])
os.makedirs(save_dir, exist_ok=True)

# Convert to a DataFrame
df = pd.DataFrame(simulated_data_neuro_musculo_state, columns=header_neuro_musculo_state)

# Save to Excel
excel_name = trial_name + '_neuro_musculo_state_simulated'
excel_path = os.path.join(save_dir, excel_name + '.xlsx')
df.to_excel(excel_path, index=False)

    # 2. Save mtu_architecture
    ##########################
# Convert to a DataFrame
df = pd.DataFrame(simulated_data_mtu_architecture, columns=header_mtu_architecture)

# Save to Excel
excel_name = trial_name + '_mtu_architecture'
excel_path = os.path.join(save_dir, excel_name + '.xlsx')
df.to_excel(excel_path, index=False)

    # 3. Save mtu_forces
    ##########################
# Convert to a DataFrame
df = pd.DataFrame(simulated_data_mtu_forces, columns=header_mtu_forces)

# Save to Excel
excel_name = trial_name + '_mtu_forces_simulated'
excel_path = os.path.join(save_dir, excel_name + '.xlsx')
df.to_excel(excel_path, index=False)
print('-------------------------------------------------------------------')
print('Process succeed')
print('simulated data saved in : ')
print(save_dir)

"""


