import os
import numpy as np
import pandas as pd
from casadi import SX, DM, vertcat, horzcat, Function, sqrt, exp, if_else, logic_and, sum1, jacobian, rootfinder, nlpsol, cos, sin, fmin, fmax, log1p,log, gradient
import math
import matplotlib.pyplot as plt
from fontTools.misc.bezierTools import epsilon
from matplotlib.widgets import Slider
import matplotlib
from live_plot_utils import OnlineCallback

matplotlib.use('TkAgg')  # or 'Qt5Agg', depending on your environment


def get_skeleton():
    """
                              # ========= 1. musculoskeletal Geometry  ========= #
    # 1.1 Topology
    B : Bones
    J : Joins
    df : degree of freedom

    Total df = 6 df
    --> J : Hip (4 df : Rotx= 0 ; Roty = 0 ; Rotz = theta_hip ; X = x ; Y = y; Z = z )
        --> B : Thigh
            --> J : knee (1 df : Rotx= 0 ; Roty = 0 ; Rotz = theta_knee ; X = 0 ; Y = 0 ; Z = 0 )
                --> B : Leg
                    --> J : Ankle (1 dF: Rotx= 0 ; Roty = 0 ; Rotz = theta_ankle ; X = 0 ; Y = 0 ; Z = 0 )
                        --> B : Talus
                            --> J : Subtalar (0 dF: Rotx= 0 ; Roty = 0 ; Rotz = 0 ; X = 0 ; Y = 0 ; Z = 0 )
                                --> B : Calc


    n_muscles = 3 ;
    NameMuscles = {'TibialisAnterior','Soleus','Gastrocnemius'} ;
    nBones = 4 ; NameBones = {'Thigh','Leg','Talus','Calcaneus'} ;
    nJoints = 4 ; NameJoints = {'Hip','Knee','Ankle','Subtalar'} ;
    """

    # ========= useful functions  ========= #
    def rototranslation_rz(translation, theta):
        # Build 4x4 homogeneous transformation matrix for rotation about Z and translation
        tx, ty, tz = translation[0], translation[1], translation[2]
        trans_matrix = SX.zeros(4, 4)
        trans_matrix[0, 0] = cos(theta)
        trans_matrix[0, 1] = -sin(theta)
        trans_matrix[1, 0] = sin(theta)
        trans_matrix[1, 1] = cos(theta)
        trans_matrix[0, 3] = tx
        trans_matrix[1, 3] = ty
        trans_matrix[2, 3] = tz
        trans_matrix[2, 2] = 1
        trans_matrix[3, 3] = 1
        return trans_matrix

    def to_global(trans_matrix, p_local):
        # Apply rotation and translation from 4x4 matrix
        return (trans_matrix @ vertcat(p_local, SX(1)))[:3]

    # Constants
    n_muscles = 3

    # ── generalised coordinates  q ∈ ℝ⁶ ───────────────────────────────────────
    x, y, z = SX.sym('x'), SX.sym('y'), SX.sym('z')
    theta_hip = SX.sym('theta_hip')
    theta_knee = SX.sym('theta_knee')
    theta_ankle = SX.sym('theta_ankle')
    q = vertcat(x, y, z, theta_hip, theta_knee, theta_ankle)

    # ── segment geometry — distal joint position in parent frame ──────────────
    thigh_geom = SX.sym('thigh_geom', 3)  # knee position in thigh frame
    leg_geom = SX.sym('leg_geom', 3)  # ankle position in leg frame
    talus_geom = SX.sym('talus_geom', 3)  # calcaneus position in talus frame
    calc_geom = SX.sym('calc_geom', 3)  # toe position in calcaneus frame
    segment_geometry = vertcat(thigh_geom, leg_geom, talus_geom, calc_geom)

    # ── muscle attachment points (local frames) ───────────────────────────────
    # tibialis anterior : leg → calcaneus, with one via point on the leg
    origin_tibialis = SX.sym('origin_tibialis', 3)
    via_tibialis = SX.sym('via_tibialis', 3)
    insertion_tibialis = SX.sym('insertion_tibialis', 3)

    # soleus : leg → calcaneus
    origin_soleus = SX.sym('origin_soleus', 3)
    insertion_soleus = SX.sym('insertion_soleus', 3)

    # gastrocnemius : thigh → calcaneus (bi-articular)
    origin_gastrocnemius = SX.sym('origin_gastrocnemius', 3)
    insertion_gastrocnemius = SX.sym('insertion_gastrocnemius', 3)

    muscle_geom = vertcat(
        origin_tibialis, origin_soleus, origin_gastrocnemius,
        insertion_tibialis, insertion_soleus, insertion_gastrocnemius,
        via_tibialis
    )

    musculoskeletal_params = vertcat(segment_geometry, muscle_geom)

    # ── forward kinematics : segment frames in global ─────────────────────────
    r_0_thigh = rototranslation_rz(vertcat(x, y, z), theta_hip)
    r_thigh_leg = rototranslation_rz(thigh_geom, -theta_knee)
    r_leg_talus = rototranslation_rz(leg_geom, theta_ankle)
    r_talus_calc = rototranslation_rz(talus_geom, 0)

    # ── express in R0 ─────────────────────────────────────────────────────────
    r_0_leg = r_0_thigh @ r_thigh_leg
    r_0_talus = r_0_leg @ r_leg_talus
    r_0_calc = r_0_talus @ r_talus_calc

    # ── muscle attachments in global frame ────────────────────────────────────
    origin_tibialis_global = to_global(r_0_leg, origin_tibialis)
    via_tibialis_global = to_global(r_0_leg, via_tibialis)
    insertion_tibialis_global = to_global(r_0_calc, insertion_tibialis)

    origin_soleus_global = to_global(r_0_leg, origin_soleus)
    insertion_soleus_global = to_global(r_0_calc, insertion_soleus)

    origin_gastrocnemius_global = to_global(r_0_thigh, origin_gastrocnemius)
    insertion_gastrocnemius_global = to_global(r_0_calc, insertion_gastrocnemius)

    # ── joint and markers ────────────────────────────────────────────────────
    origin_global = horzcat(origin_tibialis_global, origin_soleus_global, origin_gastrocnemius_global)
    insertion_global = horzcat(insertion_tibialis_global, insertion_soleus_global, insertion_gastrocnemius_global)
    via_global = horzcat(via_tibialis_global)
    hjc = r_0_thigh[0:3, 3]
    kjc = r_0_leg[0:3, 3]
    ajc = r_0_talus[0:3, 3]
    tjc = to_global(r_0_calc, calc_geom)
    calc = r_0_calc[0:3, 3]
    markers = horzcat(hjc, kjc, ajc, tjc, calc)

    musculoskeletal_states = vertcat(q, musculoskeletal_params)

    # ── MTU lengths and moment arms ───────────────────────────────────────────
    umt_length = vertcat(
        sqrt(sum1((insertion_tibialis_global - via_tibialis_global) ** 2)) +
        sqrt(sum1((via_tibialis_global - origin_tibialis_global) ** 2)),  # tibialis
        sqrt(sum1((insertion_soleus_global - origin_soleus_global) ** 2)),  # soleus
        sqrt(sum1((insertion_gastrocnemius_global - origin_gastrocnemius_global) ** 2))  # gastroc
    )

    moment_arm = jacobian(umt_length, q)

    # ── CasADi exported functions ─────────────────────────────────────────────
    forward_kinematics = Function("forward_kinematics",
                                  [musculoskeletal_states],
                                  [origin_global, insertion_global, via_global, markers],
                                  ["musculoskeletal_states"],
                                  ["origin_global", "insertion_global", "via_global", "markers"])

    get_mtu_length = Function("get_mtu_length",
                              [musculoskeletal_states],
                              [umt_length],
                              ["musculoskeletal_states"],
                              ["umt_length (tibialis, soleus, gastrocnemius)"])

    get_moment_arm = Function("get_moment_arm",
                              [musculoskeletal_states],
                              [moment_arm],
                              ["musculoskeletal_states"],
                              ["moment Arm (tibialis, soleus, gastrocnemius)"])

    return (q, moment_arm, musculoskeletal_params, musculoskeletal_states,
            forward_kinematics, get_mtu_length, get_moment_arm)

def get_fiber_active_force_length(a, normalized_fiber_length, maximal_isometric_force):
    # === Active Force-Length (S2) ===
    # First Gaussian coefficients
    b11, b21, b31, b41 = 0.814483478343008, 1.055033428970575, 0.162384573599574, 0.063303448465465
    # Second Gaussian coefficients
    b12, b22, b32, b42 = 0.433004984392647, 0.716775413397760, -0.029947116970696, 0.200356847296188
    # Third Gaussian coefficients
    b13, b23, b33, b43 = 0.100, 1.000, 0.5 * np.sqrt(0.5), 0.000

    # Assume these are defined: normalizedFiberLength, a, maximal_isometric_force
    # normalizedFiberLength should be a float or a NumPy array

    # Gaussian 1
    num1 = normalized_fiber_length - b21
    den1 = b31 + b41 * normalized_fiber_length
    fm_tilde1 = b11 * np.exp(-0.5 * (num1 ** 2) / (den1 ** 2))

    # Gaussian 2
    num2 = normalized_fiber_length - b22
    den2 = b32 + b42 * normalized_fiber_length
    fm_tilde2 = b12 * np.exp(-0.5 * (num2 ** 2) / (den2 ** 2))

    # Gaussian 3
    num3 = normalized_fiber_length - b23
    den3 = b33 + b43 * normalized_fiber_length
    fm_tilde3 = b13 * np.exp(-0.5 * (num3 ** 2) / (den3 ** 2))

    # Total normalized active force-length
    normalized_fiber_active_force_length = fm_tilde1 + fm_tilde2 + fm_tilde3

    # no forces in extrema l_f<50% = 0 [Buchanan et al. (2004) page 10
    # normalized force ∈ [0 1.7[
    normalized_fiber_active_force_length = fmax(normalized_fiber_active_force_length, 0.00)

    # Non-normalized active force
    fiber_active_force_length = a * normalized_fiber_active_force_length * maximal_isometric_force

    return fiber_active_force_length

def get_fiber_passive_force_length(normalized_fiber_length, k_fiber, maximal_isometric_force):
    # === Passive Force-Length (S3) ===
    """
        # hard max
    e0 = 0.6

    # exponetial forces (lf > lom)
    normalized_fiber_passive_force = (exp(((k_fiber * (
            normalized_fiber_length - 1)) / e0)) - 1) / (exp(k_fiber) - 1)

    # non-negative forces (lf < lom)
    normalized_fiber_passive_force = if_else(
        normalized_fiber_passive_force > 0,
        normalized_fiber_passive_force,
        0
    )

    # normalized force ∈ [0 2[
    normalized_fiber_passive_force = fmax(normalized_fiber_passive_force, 0.0)
    normalized_fiber_passive_force = fmin(normalized_fiber_passive_force, 50)

    fiber_passive_force = normalized_fiber_passive_force * maximal_isometric_force  # Non - normalized equation
    """

    # avoid inf in jackobian
    e0 = 0.6
    F_THRESHOLD = 30.0

    # ----- Point de raccord et pente sur la courbe exp -----
    arg_threshold = log(F_THRESHOLD * (exp(k_fiber) - 1.0) + 1.0)
    l_star = 1.0 + e0 / k_fiber * arg_threshold
    slope_at_threshold = (k_fiber / e0) * (F_THRESHOLD * (exp(k_fiber) - 1.0) + 1.0) / (exp(k_fiber) - 1.0)

    # ----- Branche exponentielle (arg plafonné pour éviter l'overflow,
    #       car CasADi évalue toujours les deux branches d'un if_else) -----
    arg = k_fiber * (normalized_fiber_length - 1.0) / e0
    arg_capped = fmin(arg, arg_threshold + 1.0)
    f_exp = (exp(arg_capped) - 1.0) / (exp(k_fiber) - 1.0)

    # ----- Branche logarithmique au-delà du seuil -----
    f_log = F_THRESHOLD + log(1.0 + slope_at_threshold * (normalized_fiber_length - l_star))

    # ----- Recollage C1 -----
    normalized_fiber_passive_force = if_else(
        normalized_fiber_length < l_star,
        f_exp,
        f_log
    )

    # Force positive uniquement (lf < lom)
    normalized_fiber_passive_force = fmax(normalized_fiber_passive_force, 0.0)

    fiber_passive_force = normalized_fiber_passive_force * maximal_isometric_force

    return fiber_passive_force

def get_tendon_force_length(normalized_tendon_length, k_tendon, maximal_isometric_force):
    # Tendon force-length (S1)
    c1 = 0.200;     c2 = 0.995;     c3 = 0.250  # tendon parameters

    # --- Calcul de f(1) pour centrer la courbe ---
    arg_at_1 = k_tendon * (1.0 - c2)
    arg_safe_at_1 = 50.0 - log1p(exp(50.0 - arg_at_1))
    f_at_1 = c1 * exp(arg_safe_at_1) - c3

    # --- Courbe principale ---
    arg = k_tendon * (normalized_tendon_length - c2)
    arg_safe = 50.0 - log1p(exp(50.0 - arg))
    normalized_tendon_force_curve = c1 * exp(arg_safe) - c3

    # --- Condition f(1) = 0 ---
    normalized_tendon_force_curve = normalized_tendon_force_curve - f_at_1

    normalized_tendon_force_curve = fmax(normalized_tendon_force_curve, 0)
    tendon_force = (normalized_tendon_force_curve * maximal_isometric_force) # /k_tendon

    return tendon_force

def get_muscle_total_force(fiber_active_force_length, normalized_fiber_force_velocity, fiber_passive_force):
    return fiber_active_force_length * normalized_fiber_force_velocity + fiber_passive_force

def get_muscle_dynamic(q, moment_arm, musculoskeletal, n_muscles,
                       param_config=None, fixed_values=None):
    """
    Build CasADi functions for muscle-tendon dynamics (De Groote 2016 model),
    including equilibrium rootfinders (per-muscle and all-muscles).

    Muscle-tendon parameters (ℓom, φo, Fom, ℓst, kt, km) can be declared either
    as symbolic variables (for optimization / identification) or hardcoded as
    numerical constants, independently for each parameter.

    Parameters
    ----------
    q, musculoskeletal : CasADi SX
    n_muscles : int
    param_config : dict {param_name: 'sym' | 'fixed'}, optional
        Default : all 'sym'.
    fixed_values : dict {param_name: scalar or array-like (n_muscles,)}, optional
        Required for every parameter marked 'fixed'.

    Returns
    -------
    dict with CasADi Functions and rootfinders.

    References
    ----------
    De Groote F, Kinney AL, Rao A V, Fregly BJ. Ann Biomed Eng. 2016;
    44: 2922–2936. doi:10.1007/s10439-016-1591-9
    """

    # ========= 1. Résolution de la configuration des paramètres ========= #
    MTU_PARAM_SPEC = [
        ('l0m', 'ℓom'),
        ('phi0', 'φo'),
        ('f0m', 'Fom'),
        ('km', 'km'),
        ('lst', 'ℓst'),
        ('kt', 'kt'),
    ]
    valid_keys = [k for k, _ in MTU_PARAM_SPEC]

    if param_config is None:
        param_config = {k: 'sym' for k in valid_keys}
    fixed_values = fixed_values or {}

    for k in valid_keys:
        if k not in param_config:
            raise ValueError(
                f"param_config missing entry for '{k}'. Required keys: {valid_keys}"
            )
        if param_config[k] not in ('sym', 'fixed'):
            raise ValueError(
                f"param_config['{k}'] must be 'sym' or 'fixed', got {param_config[k]!r}"
            )

    for k, status in param_config.items():
        if status == 'fixed' and k not in fixed_values:
            raise ValueError(
                f"Parameter '{k}' is marked 'fixed' but missing from fixed_values."
            )

    # ========= 2. Construction des variables MTU ========= #
    def _build_mtu_param(key, display_name):
        if param_config[key] == 'sym':
            return SX.sym(display_name, n_muscles)
        val = np.asarray(fixed_values[key], dtype=float).reshape(-1)
        if val.size == 1:
            val = np.full(n_muscles, val.item())
        if val.size != n_muscles:
            raise ValueError(
                f"fixed_values['{key}'] has size {val.size}, expected 1 or {n_muscles}."
            )
        return SX(val)

    mtu_vars = {key: _build_mtu_param(key, name) for key, name in MTU_PARAM_SPEC}

    optimal_fiber_length = mtu_vars['l0m']
    phi0 = mtu_vars['phi0']
    maximal_isometric_force = mtu_vars['f0m']
    tendon_slack_length = mtu_vars['lst']
    k_tendon = mtu_vars['kt']
    k_fiber = mtu_vars['km']

    # Vecteur des paramètres SYMBOLIQUES (pour optim / rootfinder)
    sym_keys = [k for k, _ in MTU_PARAM_SPEC if param_config[k] == 'sym']
    sym_params = [mtu_vars[k] for k in sym_keys]
    if sym_params:
        muscle_tendon_parameters = vertcat(*sym_params)
    else:
        muscle_tendon_parameters = SX.sym('mtu_params_empty', 0)

    # Version "single muscle" : premier élément de chaque paramètre symbolique
    sym_params_single = [mtu_vars[k][0] for k in sym_keys]
    if sym_params_single:
        muscle_tendon_parameters_single = vertcat(*sym_params_single)
    else:
        muscle_tendon_parameters_single = SX.sym('mtu_params_single_empty', 0)

    # ========= 3. États muscle-tendon ========= #
    a = SX.sym('muscle_activation', n_muscles)
    fiber_length = SX.sym('fiber_length', n_muscles)
    tendon_length = SX.sym('tendon_length', n_muscles)
    pennation_angle = SX.sym('pennation_angle', n_muscles)

    rooted_variables = vertcat(fiber_length, pennation_angle, tendon_length)

    # ========= 4. Équations MTU (De Groote 2016) ========= #
    normalized_tendon_length = tendon_length / tendon_slack_length
    normalized_fiber_length = fiber_length / optimal_fiber_length

    tendon_force = get_tendon_force_length(
        normalized_tendon_length, k_tendon, maximal_isometric_force
    )
    fiber_active_force_length = get_fiber_active_force_length(
        a, normalized_fiber_length, maximal_isometric_force
    )
    fiber_passive_force = get_fiber_passive_force_length(
        normalized_fiber_length, k_fiber, maximal_isometric_force
    )
    normalized_fiber_force_velocity = 1.0  # quasi-static
    fiber_force = get_muscle_total_force(
        fiber_active_force_length,
        normalized_fiber_force_velocity,
        fiber_passive_force,
    )

    # ========= 5. Wrap into CasADi functions ========= #
    neuromusculoskeletal_state = vertcat(a, q, musculoskeletal)
    all_states = vertcat(neuromusculoskeletal_state, rooted_variables)

    # ========= 6. Rootfinder : équilibre muscle-tendon ========= #
    # Input géométrique : longueur MTU (ℓmt) — donnée par la cinématique
    l_mtu = SX.sym('umt_length', n_muscles)

    # --- Contraintes d'équilibre (De Groote 2016, eqs S5-S7) --- #
    # g5 : ℓmt = cos(φ) · ℓm + ℓt          (longueur MTU)
    # g6 : ℓom · sin(φo) = ℓm · sin(φ)     (conservation du volume musculaire)
    # g7 : Ff · cos(φ) = Ft                (équilibre des forces axiales)
    """
    g5 = l_mtu - (cos(pennation_angle) * fiber_length + tendon_length)
    g6 = (optimal_fiber_length * sin(phi0)) - (fiber_length * sin(pennation_angle))
    g7 = fiber_force * cos(pennation_angle) - tendon_force
    """

    g5 = (l_mtu - (cos(pennation_angle) * fiber_length + tendon_length)) / optimal_fiber_length
    g6 = (optimal_fiber_length * sin(phi0) - fiber_length * sin(pennation_angle)) / optimal_fiber_length
    g7 = (fiber_force * cos(pennation_angle) - tendon_force) / maximal_isometric_force


    # --- Problème single-muscle : 3 inconnues, 3 équations --- #
    unknown_single = vertcat(fiber_length[0], pennation_angle[0], tendon_length[0])
    known_single = vertcat(a[0], l_mtu[0], muscle_tendon_parameters_single)

    equilibrium_error_single_muscle = Function(
        'equilibrium_error_single_muscle',
        [unknown_single, known_single],
        [vertcat(g5[0], g6[0], g7[0])],
        ['x', 'p'],
        ['residuals'],
    )


    opts_newton_single = {
        "abstol": 1e-8,
        "max_iter": 5000,
        "error_on_fail": False,
        "print_iteration": False,
    }

    equilibrate_muscle_tendon_single_muscle = rootfinder(
        'equilibrate_muscle_tendon_single_muscle',
        'newton',
        equilibrium_error_single_muscle,
        opts_newton_single,
    )

    # --- Problème all-muscles : 3·n_muscles inconnues, 3·n_muscles équations --- #
    unknown_all = vertcat(fiber_length, pennation_angle, tendon_length)
    known_all = vertcat(a, l_mtu, muscle_tendon_parameters)

    equilibrium_error_all_muscle = Function(
        'equilibrium_error_all_muscle',
        [unknown_all, known_all],
        [vertcat(g5, g6, g7)],
        ['x', 'p'],
        ['residuals'],
    )

    opts_newton_all = {
        "abstol": 1e-8,
        "max_iter": 5000,
        "error_on_fail": False,
        "print_iteration": False,
    }
    equilibrate_muscle_tendon_all = rootfinder(
        'equilibrate_muscle_tendon_all',
        'newton',
        equilibrium_error_all_muscle,
        opts_newton_all,
    )

    # ========= 4. Computing Joint Moments and Angles    ========= #
    joint_torque = moment_arm * tendon_force
    joint_torque = sum1(joint_torque[:, -1:])

    get_tendon_force_from_tendon_length = Function(
        'get_tendon_force_from_tendon_length',
        [tendon_length, muscle_tendon_parameters],
        [tendon_force],
        ['tendon_length', 'muscle_tendon_parameters'],
        ['tendon_force']
    )

    get_fiber_force_from_fiber_length = Function(
        'get_fiber_force_from_fiber_length',
        [a, fiber_length, muscle_tendon_parameters],
        [fiber_force],
        ['muscle_activation', 'fiber_length', 'muscle_tendon_parameters'],
        ['fiber_force']
    )

    def _wrap(name, output, output_name):
        return Function(
            name,
            [all_states, muscle_tendon_parameters],
            [output],
            ['all_states', 'muscle_tendon_parameters'],
            [output_name],
        )

    get_tendon_force_fn = _wrap('get_tendon_force', tendon_force, 'tendon_force')
    get_fiber_force_fn = _wrap('get_fiber_force', fiber_force, 'fiber_force')
    get_fiber_active_force_fn = _wrap('get_fiber_active_force', fiber_active_force_length, 'fiber_active_force')
    get_fiber_passive_force_fn = _wrap('get_fiber_passive_force', fiber_passive_force, 'fiber_passive_force')
    get_joint_moment_fn = _wrap('get_joint_moment', joint_torque, 'joint_torque')

    # ========= 7. Assemblage du retour ========= #
    return {
        # Fonctions representation
        'get_tendon_force_from_tendon_length': get_tendon_force_from_tendon_length,
        'get_fiber_force_from_fiber_length': get_fiber_force_from_fiber_length,
        # Fonctions de force
        'get_tendon_force': get_tendon_force_fn,
        'get_fiber_force': get_fiber_force_fn,
        'get_fiber_active_force': get_fiber_active_force_fn,
        'get_fiber_passive_force': get_fiber_passive_force_fn,
        'get_joint_moment': get_joint_moment_fn,
        # Résidus d'équilibre (utiles pour diagnostic)
        'equilibrium_error_single_muscle': equilibrium_error_single_muscle,
        'equilibrium_error_all_muscle': equilibrium_error_all_muscle,
        # Rootfinders
        'equilibrate_muscle_tendon_single_muscle': equilibrate_muscle_tendon_single_muscle,
        'equilibrate_muscle_tendon_all': equilibrate_muscle_tendon_all,
        # Métadonnées
        'muscle_tendon_parameters': muscle_tendon_parameters,
        'muscle_tendon_parameters_single': muscle_tendon_parameters_single,
        'sym_param_order': sym_keys,
        'param_config': dict(param_config),
    }

def get_model_equation(param_config=None, fixed_values=None, n_muscles=3):
    """
    Build the full musculoskeletal model by merging skeleton kinematics
    and muscle-tendon dynamics into a unified set of CasADi functions.

    Parameters
    ----------
    param_config : dict, optional
        Per-parameter status ('sym' or 'fixed') for the 6 MTU parameters
        (l0m, phi0, f0m, lst, kt, km). Default : geometric params symbolic
        (l0m, phi0, f0m, lst), shape params fixed (kt, km).
    fixed_values : dict, optional
        Numerical values for parameters marked 'fixed'. Scalars are broadcast
        to all n_muscles. Default : kt=3.0, km=4.0.
    n_muscles : int, optional
        Number of muscles in the model. Default: 3 (tib_ant, soleus, gast).

    Returns
    -------
    casadi_functions : dict[str, casadi.Function]
        Unified collection of all CasADi functions describing the model :
        kinematics + muscle-tendon dynamics + equilibrium rootfinders.
    sym_param_order : list of str
        Ordered list of MTU parameter names that are symbolic, indicating
        the layout of the 'muscle_tendon_parameters' input vector.
    definition : dict
        Documentation of the symbolic variables used in the model.
    """

    # ========= 1. Default settings ========= #
    if param_config is None:
        param_config = {
            'l0m': 'sym',
            'phi0': 'sym',
            'f0m': 'sym',
            'km': 'fixed',
            'lst': 'sym',
            'kt': 'fixed',
        }
    if fixed_values is None:
        fixed_values = {
            'kt': 35,
            'km': 4.0,
        }

    # ========= 2. Skeleton (kinematics) ========= #
    (q, moment_arm, musculoskeletal_params, musculoskeletal_states,
     forward_kinematics, get_mtu_length, get_moment_arm) = get_skeleton()

    # ========= 3. Muscle-tendon dynamics ========= #
    funcs = get_muscle_dynamic(
        q, moment_arm, musculoskeletal_params,
        n_muscles=n_muscles,
        param_config=param_config,
        fixed_values=fixed_values,
    )

    # ========= 4. Merge all CasADi functions ========= #
    casadi_functions = {
        name: obj for name, obj in funcs.items()
        if isinstance(obj, Function)
    }

    extras = {
        'forward_kinematics': forward_kinematics,
        'get_mtu_length': get_mtu_length,
        'get_moment_arm': get_moment_arm,
        'mtu_parameters_sym': funcs['muscle_tendon_parameters'],
        'mtu_parameters_name': funcs['sym_param_order'],
    }

    name_clash = set(casadi_functions) & set(extras)
    if name_clash:
        raise ValueError(f"Key collision in casadi_functions: {name_clash}")
    casadi_functions.update(extras)

    # ========= 5. Documentation des variables symboliques ========= #
    definition = {
        'a': "Neuromuscular activation (between 0 and 1)",
        'q': "Spatial skeleton configuration "
             "(x, y, z, theta_hip, theta_knee, theta_ankle)",
        'musculoskeletal': "Musculoskeletal parameters "
                           "(segment_geometry, muscle_insertion, "
                           "Local_ViaPoint_tibialis_anterior)",
        'neuromusculoskeletal_state': "Concatenation (a, q, musculoskeletal)",
        'muscle_tendon_parameters': f"Symbolic MTU parameters in order: "
                                    f"{funcs['sym_param_order']}",
        'rooted_variables': "Hidden states resolved by rootfinder "
                            "(fiber_length, pennation_angle, tendon_length)",
        'l_mtu': "Total muscle-tendon unit length (input to equilibrium)",
        'p': "Full input vector for force functions: "
             "(neuromusculoskeletal_state, muscle_tendon_parameters)",
    }

    return casadi_functions, funcs['muscle_tendon_parameters'], definition

def test_model(skeleton_num, muscle_tendon_parameters_num, casadi_function, a_num, q_num):
    # === plot musculosquelton system ===
    origin_global_num, insertion_global_num, via_global, markers_num = casadi_function['forward_kinematics'](
        np.concatenate((q_num, skeleton_num)).reshape(-1, 1))

    # === get muscle length ===
    mtu_length = casadi_function['get_mtu_length'](np.concatenate((q_num, skeleton_num)).reshape(-1, 1))

    tibialis_length = float(mtu_length[0])
    soleus_length = float(mtu_length[1])
    gastrocnemius_length = float(mtu_length[2])

    print('\n', '================== mtu length length ==================')
    print('tibialis: ', tibialis_length, 'm')
    print('soleus: ', soleus_length, 'm')
    print('gastrocnemius: ', gastrocnemius_length, 'm')

    musculoskeletal_states_num = np.concatenate((q_num, skeleton_num))
    neuromusculoskeletal_state_num = np.concatenate((a_num, musculoskeletal_states_num))

    # === get muscle moment arm ===
    mtu_moment_arm = casadi_function['get_moment_arm'](musculoskeletal_states_num)

    tibialis_moment_arm = float(mtu_moment_arm[0, -1])
    soleus_moment_arm = float(mtu_moment_arm[1, -1])
    gastrocnemius_moment_arm = float(mtu_moment_arm[2, -1])

    print('\n', '================== mtu moment arm ==================')
    print('tibialis: ', tibialis_moment_arm, 'm-1')
    print('soleus: ', soleus_moment_arm, 'm-1')
    print('gastrocnemius: ', gastrocnemius_moment_arm, 'm-1')

    # === muscle equilibrium ===

    """
    muscle_tendon_parameters_num = np.array(muscle_tendon_parameters_num)

    n_muscle = 3
    n_param = len(muscle_tendon_parameters_num)/n_muscle
    muscle_tendon_parameters_num_tibialis = muscle_tendon_parameters_num[[0, 3, 6, 9]]
    x_opt_tibialis,state_tibialis = root_muscle_dynamics(a_num[0], float(mtu_length[0]),muscle_tendon_parameters_num[[0, 3, 6, 9]], 'tibialis', casadi_function)
    x_opt_soleus,state_soleus = root_muscle_dynamics(a_num[1], float(mtu_length[1]), muscle_tendon_parameters_num[[1, 4, 7, 10]], 'soleus', casadi_function)
    x_opt_gastrocnemius,state_gastrocnemius = root_muscle_dynamics(a_num[2], float(mtu_length[2]),muscle_tendon_parameters_num[[3, 5, 8, 11]], 'gastrocnemius',casadi_function)

    # rooted_variables = [fiber length, pennation angle and tendon length]
    rooted_variables = np.array([x_opt_tibialis[0, 0], x_opt_soleus[0, 0], x_opt_gastrocnemius[0, 0],
                                 x_opt_tibialis[1, 0], x_opt_soleus[1, 0], x_opt_gastrocnemius[1, 0],
                                 x_opt_tibialis[2, 0], x_opt_soleus[2, 0], x_opt_gastrocnemius[2, 0]]).flatten()
    """

    results = solve_all_muscles(a_num, mtu_length, muscle_tendon_parameters_num, casadi_function)

    x_opt_tibialis = results['tibialis']['x_opt']
    x_opt_soleus = results['soleus']['x_opt']
    x_opt_gastrocnemius = results['gastrocnemius']['x_opt']
    # rooted_variables = [fiber length, pennation angle and tendon length]

    rooted_variables = np.array([x_opt_tibialis[0, 0], x_opt_soleus[0, 0], x_opt_gastrocnemius[0, 0],
                                 x_opt_tibialis[1, 0], x_opt_soleus[1, 0], x_opt_gastrocnemius[1, 0],
                                 x_opt_tibialis[2, 0], x_opt_soleus[2, 0], x_opt_gastrocnemius[2, 0]]).flatten()

    all_state = np.concatenate([neuromusculoskeletal_state_num, rooted_variables])

    tendon_force = casadi_function['get_tendon_force'](all_state, muscle_tendon_parameters_num)
    fiber_force = casadi_function['get_fiber_force'](all_state, muscle_tendon_parameters_num)
    fiber_passive_force = casadi_function['get_fiber_passive_force'](all_state, muscle_tendon_parameters_num)
    fiber_active_force = casadi_function['get_fiber_active_force'](all_state, muscle_tendon_parameters_num)

    ankle_torque = casadi_function['get_joint_moment'](all_state, muscle_tendon_parameters_num)

    print('\n================== simulatied ==================')
    print('\n   mtu architecture:')
    print('     - fiber length:')
    print(f'tibialis: {float(x_opt_tibialis[0, 0]):.4f} m')
    print(f'soleus: {float(x_opt_soleus[0, 0]):.4f} m')
    print(f'gastrocnemius: {float(x_opt_gastrocnemius[0, 0]):.4f} m')
    print('     - pennation angle:')
    print(f'tibialis: {float(np.rad2deg(float(np.array(x_opt_tibialis[1, 0]).flatten()[0]))):.4f} deg')
    print(f'soleus: {float(np.rad2deg(float(np.array(x_opt_soleus[1, 0]).flatten()[0]))):.4f} deg')
    print(f'gastrocnemius: {float(np.rad2deg(float(np.array(x_opt_gastrocnemius[1, 0]).flatten()[0]))):.4f} deg')
    print('     - tendon length:')
    print(f'tibialis: {float(x_opt_tibialis[2, 0]):.4f} m')
    print(f'soleus: {float(x_opt_soleus[2, 0]):.4f} m')
    print(f'gastrocnemius: {float(x_opt_gastrocnemius[2, 0]):.4f} m')

    print('\n   tendon force:')
    print(f'tibialis: {float(tendon_force[0]):.4f} N')
    print(f'soleus: {float(tendon_force[1]):.4f} N')
    print(f'gastrocnemius: {float(tendon_force[2]):.4f} N')

    print('\n   fiber force:')
    print('     - fiber force:')
    print(f'tibialis: {float(fiber_force[0]):.4f} N')
    print(f'soleus: {float(fiber_force[1]):.4f} N')
    print(f'gastrocnemius: {float(fiber_force[2]):.4f} N')

    print('     - fiber passive force:')
    print(f'tibialis: {float(fiber_passive_force[0]):.4f} N')
    print(f'soleus: {float(fiber_passive_force[1]):.4f} N')
    print(f'gastrocnemius: {float(fiber_passive_force[2]):.4f} N')

    print('     - fiber active force:')
    print(f'tibialis: {float(fiber_active_force[0]):.4f} N')
    print(f'soleus: {float(fiber_active_force[1]):.4f} N')
    print(f'gastrocnemius: {float(fiber_active_force[2]):.4f} N')

    print('\n   joint torque:')
    print(f'ankle moment: {float(ankle_torque[0]):.4f} N.m')

def solve_all_muscles(a_num, mtu_length, mtu_params, casadi_function):
    """
    Solve the muscle-tendon equilibrium for every muscle.

    mtu_params is assumed to be in PARAM-MAJOR layout :
        [p1_m1, p1_m2, ..., p2_m1, p2_m2, ...]
    For muscle index m, parameters are at indices m, m+n_muscles, m+2*n_muscles...

    Returns
    -------
    dict {muscle_name: {'x_opt': ..., 'state': ..., 'params': ...}}
    """
    muscle_names = 'tibialis', 'soleus', 'gastrocnemius'
    mtu_params = np.asarray(mtu_params)
    n_muscles = len(muscle_names)

    results = {}
    for m_idx, m_name in enumerate(muscle_names):
        params_m = mtu_params[m_idx::n_muscles]
        x_opt, state, residuals = root_muscle_dynamics(
            a_num[m_idx],
            float(mtu_length[m_idx]),
            params_m,
            m_name,
            casadi_function,
        )
        results[m_name] = {'x_opt': x_opt, 'state': state, 'params': params_m,'residuals':residuals}
    return results

def root_muscle_dynamics(a, lmtu, parameters, muscle_name, casadi_function):
    """
    Solves for the muscle-tendon equilibrium state using root-finding with CasADi.

    Parameters:
    - a (float): Muscle activation level.
    - lmtu (float): Muscle-tendon unit length.
    - parameters (np.ndarray): Array of muscle-specific parameters required by the model.
    - muscle_name (str): Name of the muscle (used for print/logging purposes).
    - casadi_function (dict): Dictionary containing CasADi functions:
        - 'equilibrate_muscle_tendon_single_muscle': Root-finding function to solve muscle-tendon equilibrium.
        - 'equilibrium_error_single_muscle': Function to evaluate equilibrium error (residuals).

    Returns:
    - x_opt (np.ndarray): Optimized muscle-tendon state vector satisfying equilibrium conditions.
    - equilibrium_status (str): Status indicating whether equilibrium was successfully achieved ('success' or 'fail_best_effort').
    - residuals_print: Residuals associated with the returned x_opt.

    Description:
    This function attempts to solve the equilibrium condition of a single muscle-tendon unit
    using a root-finding algorithm. It allows up to `max_attempts` (default: 1000) to find a
    solution where all residuals fall below a defined tolerance (`lim_residuals`).

    A "memory" mechanism keeps track of the best attempt so far (i.e., the one with the
    smallest maximum residual). If no attempt satisfies the tolerance within `max_attempts`,
    the best solution found is returned with status 'fail_best_effort'.
    """
    lim_residuals = 1e-9
    equilibrium_status = 'fail'
    max_attempts = 50
    attempt = 0

    parameters_name = casadi_function['mtu_parameters_name']
    required = 'l0m', 'phi0', 'f0m', 'lst'  # extract (ℓom, φo, Fom, ℓst)
    params = get_named_params(parameters, parameters_name, required)

    # --- Memory: best attempt so far ---
    best_x_opt = None
    best_residuals_print = None
    best_residuals_metric = np.inf  # metric used for comparison (max abs residual)

    while equilibrium_status == 'fail' and attempt < max_attempts:
        attempt += 1
        x_start = x_start_equi_tendon(a, lmtu, params)
        print(x_start)

        # Reset locals to detect failures cleanly
        x_opt = None
        residuals = None
        residuals_print = None

        try:
            x_opt = casadi_function['equilibrate_muscle_tendon_single_muscle'](
                x_start,
                np.array([a, lmtu] + parameters.tolist())
            )

            if x_opt[1] < 0 or x_opt[1] > math.pi / 2:
                x_opt[1] = DM(float(x_opt[1]) % math.pi / 2)
                x_opt[0] = DM(abs(float(x_opt[0])))
                x_opt[2] = lmtu - (np.cos(x_opt[1]) * x_opt[0])

            if x_opt[2] < params[3]: # tendon length < tsl so tendon length = tsl
                l0m = float(params[0])
                phi0 = float(params[1])
                lst = float(params[3])

                # constante de volume (invariant g6) : produit l_f · sin(φ) constant
                # à l'état de référence (l0m, phi0)
                h_muscle = l0m * np.sin(phi0)  # "hauteur" anatomique du muscle

                tendon_length = lst
                l_f_along_mtu = lmtu - tendon_length
                pennation_angle = np.arctan2(h_muscle, l_f_along_mtu)
                fiber_length = l_f_along_mtu / np.cos(pennation_angle)

                x_opt[0] = fiber_length
                x_opt[1] = pennation_angle
                x_opt[2] = tendon_length


            residuals_print = casadi_function['equilibrium_error_single_muscle'](
                x_opt,
                np.array([a, lmtu] + parameters.tolist())
            )

            residuals = np.array(residuals_print).flatten()

        except Exception as e:
            print("Something went wrong:", e)
            continue  # skip memory update if the call itself failed

        # --- Update memory with the best attempt so far ---
        # Use max absolute residual as comparison metric; ignore NaNs
        if residuals is not None and not np.any(np.isnan(residuals)):
            current_metric = np.max(np.abs(residuals))
            if current_metric < best_residuals_metric:
                best_residuals_metric = current_metric
                best_x_opt = x_opt
                best_residuals_print = residuals_print

        if np.any((residuals > lim_residuals) | np.isnan(residuals)):
            print(f"[Attempt {attempt}] Root-finding failed: residual too high ({residuals})")
        else:
            print(f"[Attempt {attempt}] Success: residual within tolerance")
            equilibrium_status = 'success'

    # --- Fallback: if no success, return the best attempt found ---
    if equilibrium_status == 'fail':
        if best_x_opt is not None:
            x_opt = best_x_opt
            residuals_print = best_residuals_print
            print(f"\n[!] No attempt reached tolerance ({lim_residuals}).")
            print(f"    Returning best attempt with max|residual| = {best_residuals_metric:.3e}")
            equilibrium_status = 'success'
        else:
            print(f"\n[!] No valid attempt found within {max_attempts} tries.")


    print('\n', '=============== dynamics: ', muscle_name, ' ===============')
    print('residuals: length equilibrium, architecture equilibrium and force equilibrium')
    print('residuals ', residuals_print)
    print('x_opt: fiber length, pennation angle and tendon lenght')
    print('x_opt ', x_opt)
    print('equilibrium status: ', equilibrium_status)
    print('=============== ========================= ===============', '\n', '\n', '\n', )
    return x_opt, equilibrium_status, residuals_print

def x_start_equi_tendon(a, lmtu, parameters, rng=None):
    """
    Initialisation SIMPLE et géométrique de l'équilibre musculo-tendineux.

    Principe : on impose le tendon, puis on déduit (fibre, pennation) pour
    satisfaire EXACTEMENT g5 et g6. g7 (équilibre des forces) sera résolu
    ensuite par le rootfinder.

        ℓt  : imposé      = max(1.03·ℓst, ℓst)   (tendon légèrement étiré,
                                                   jamais sous sa longueur repos)
        g5  : ℓf·cos(φ) = ℓmtu − ℓt              (projection axiale)
        g6  : ℓf·sin(φ) = ℓom·sin(φo) = h        (hauteur constante)

        ⇒ triangle rectangle :
            φ  = arctan( h / (ℓmtu − ℓt) )
            ℓf = sqrt( (ℓmtu − ℓt)² + h² )

    Parameters
    ----------
    a : float            -- activation ∈ [0, 1] (non utilisée ici, init géométrique)
    lmtu : float         -- longueur courante du MTU
    parameters : (4,)    -- (ℓom, φo, F0m, ℓst)
    rng : ignoré (gardé pour compatibilité de signature)

    Returns
    -------
    x_start : np.ndarray (3,) -- [fiber_length, pennation_angle, tendon_length]
    """
    a    = float(a)
    l0m  = float(parameters[0])
    phi0 = float(parameters[1])
    lst  = float(parameters[3])

    # --- Tendon piloté par l'activation : repos = ℓst, +3% d'étirement à a=1 ---
    # ⚠ garde 0.03 (physiologique, ℓt/ℓst ∈ [1.0, 1.03]).
    #   a * 1.03 donnerait jusqu'à 103% d'étirement → overflow de l'exp tendon.
    tendon_length = lst + a * 0.03 * lst

    # --- Hauteur du muscle (invariant g6) ---
    h_muscle = l0m * np.sin(phi0)

    # --- Projection axiale de la fibre (depuis g5) ---
    l_f_axial = lmtu - tendon_length

    # Sécurité : si le tendon "mange" tout le MTU, la projection devient
    # négative/nulle → repli sur une fibre courte cohérente.
    if l_f_axial <= 1e-6:
        l_f_axial = 0.5 * l0m

    # --- Pennation et longueur de fibre (triangle rectangle g5+g6) ---
    pennation_angle = np.arctan2(h_muscle, l_f_axial)
    fiber_length = np.sqrt(l_f_axial**2 + h_muscle**2)

    return np.array([fiber_length, pennation_angle, tendon_length])


def x_start_equi_fiber(a, lmtu, parameters, rng=None):
    """
    Consistent initialization for the muscle-tendon equilibrium constraints.

    "Fiber" variant: starts from muscle fiber shortening with activation,
    then derives the pennation angle (g6) and tendon length (g5) by
    construction.

    Guarantees:
        - fiber_length > 0  and  l_f >= h_muscle (consistency with g6)
        - 0 < pennation_angle < pi/2
        - tendon_length > 0
        - g5 (MTU length) satisfied by construction
        - g6 (constant volume) satisfied by construction

    Equilibrium constraints:
        g5: LMTU = cos(phi) * l_f + l_t
        g6: l0m * sin(phi0) = l_f * sin(phi)
        g7: F_M * cos(phi) = F_T

    Muscle-Tendon Parameters: (l0m, phi0, F0m, lst)

    Parameters
    ----------
    a : float
        Current muscle activation (used to shorten the fiber).
    lmtu : float
        Current muscle-tendon unit length.
    parameters : array-like, shape (4,)
        (l0m, phi0, F0m, lst).
    rng : np.random.Generator, optional
        For reproducibility.

    Returns
    -------
    x_start : np.ndarray, shape (3,)
        [fiber_length, pennation_angle, tendon_length]
    """
    rng = rng if rng is not None else np.random.default_rng()

    l0m = float(parameters[0])
    phi0 = float(parameters[1])
    lst = float(parameters[3])

    # ----- Bornes physiques sur l'angle de pennation -----
    PHI_MIN = np.deg2rad(1.0)  # évite tan(0) et division par zéro
    PHI_MAX = np.deg2rad(85.0)  # évite cos(phi) -> 0

    # ----- Invariant de volume (g6) : l_f * sin(phi) = cste -----
    h_muscle = l0m * np.sin(phi0)  # "hauteur" anatomique du muscle

    # Bornes admissibles sur l_f compatibles avec g6 et [PHI_MIN, PHI_MAX]
    l_f_min = h_muscle / np.sin(PHI_MAX)
    l_f_max = h_muscle / np.sin(PHI_MIN)

    # ----- Initial guess géométrique tenant compte de a et de l_mt -----
    # Hypothèse : déformation tendineuse croissante avec l'activation
    # (tendon ~rigide au repos, ~4% d'étirement à activation maximale)
    epsilon_t_max = 0.03
    l_t_guess = lst * (1.0 + epsilon_t_max * a)

    # Longueur de fibre projetée sur l'axe du tendon (clippée pour éviter
    # les valeurs négatives en cas de MTU très court)
    l_f_along = np.maximum(lmtu - l_t_guess, 0.0)

    # Longueur de fibre totale via g6 : l_f^2 = (l_mt - l_t)^2 + h^2
    fiber_length = np.sqrt(l_f_along ** 2 + h_muscle ** 2)

    # Petite perturbation aléatoire pour éviter de démarrer toujours
    # au même point (utile si tu fais du multi-start)
    fiber_length *= (1.0 + rng.uniform(-0.15, 0.15))

    # Sécurité : on reste dans les bornes admissibles de g6
    fiber_length = np.clip(fiber_length, l_f_min, l_f_max)

    # 2. derive phi from g6: sin(phi) * l_f = h_muscle
    sin_phi = h_muscle / fiber_length
    pennation_angle = np.arcsin(sin_phi)
    pennation_angle = np.clip(pennation_angle, PHI_MIN, PHI_MAX)

    # 3. derive l_t from g5: l_t = LMTU - cos(phi) * l_f
    tendon_length = lmtu - np.cos(pennation_angle) * fiber_length

    return np.array([fiber_length, pennation_angle, tendon_length])

def x_start_equi_mixte(a, lmtu, parameters, rng=None):
    """
    Consistent initialization for the muscle-tendon equilibrium constraints.

    "Fiber" variant: starts from muscle fiber shortening with activation,
    then derives the pennation angle (g6) and tendon length (g5) by
    construction.

    Guarantees:
        - fiber_length > 0  and  l_f >= h_muscle (consistency with g6)
        - 0 < pennation_angle < pi/2
        - tendon_length > 0
        - g5 (MTU length) satisfied by construction
        - g6 (constant volume) satisfied by construction

    Equilibrium constraints:
        g5: LMTU = cos(phi) * l_f + l_t
        g6: l0m * sin(phi0) = l_f * sin(phi)
        g7: F_M * cos(phi) = F_T

    Muscle-Tendon Parameters: (l0m, phi0, F0m, lst)

    Parameters
    ----------
    a : float
        Current muscle activation (used to shorten the fiber).
    lmtu : float
        Current muscle-tendon unit length.
    parameters : array-like, shape (4,)
        (l0m, phi0, F0m, lst).
    rng : np.random.Generator, optional
        For reproducibility.

    Returns
    -------
    x_start : np.ndarray, shape (3,)
        [fiber_length, pennation_angle, tendon_length]
    """
    rng = rng if rng is not None else np.random.default_rng()

    l0m = float(parameters[0])
    phi0 = float(parameters[1])
    lst = float(parameters[3])

    # ----- Bornes physiques sur l'angle de pennation -----
    PHI_MIN = np.deg2rad(1.0)  # évite tan(0) et division par zéro
    PHI_MAX = np.deg2rad(85.0)  # évite cos(phi) -> 0

    # ----- Invariant de volume (g6) : l_f * sin(phi) = cste -----
    h_muscle = l0m * np.sin(phi0)  # "hauteur" anatomique du muscle

    # Bornes admissibles sur l_f compatibles avec g6 et [PHI_MIN, PHI_MAX]
    l_f_min = h_muscle / np.sin(PHI_MAX)
    l_f_max = h_muscle / np.sin(PHI_MIN)

    # ----- is the lumt stretched -----
    slack_lmtu = np.cos(phi0) * l0m + lst

    if slack_lmtu < lmtu: # is  stretched
        # guess by the tendon lenthening
        # ----- Initial guess géométrique tenant compte de a et de l_mt -----
        # Hypothèse : déformation tendineuse croissante avec l'activation
        # (tendon ~rigide au repos, ~4% d'étirement à activation maximale)
        epsilon_t_max = 0.03
        l_t_guess = lst * (1.0 + epsilon_t_max * a)

        # Longueur de fibre projetée sur l'axe du tendon (clippée pour éviter
        # les valeurs négatives en cas de MTU très court)
        l_f_along = np.maximum(lmtu - l_t_guess, 0.0)

        # Longueur de fibre totale via g6 : l_f^2 = (l_mt - l_t)^2 + h^2
        fiber_length = np.sqrt(l_f_along ** 2 + h_muscle ** 2)

        # Petite perturbation aléatoire pour éviter de démarrer toujours
        # au même point (utile si tu fais du multi-start) 5%
        fiber_length *= (1.0 + rng.uniform(-0.005, 0.005))

        # Sécurité : on reste dans les bornes admissibles de g6
        fiber_length = np.clip(fiber_length, l_f_min, l_f_max)

        # 2. derive phi from g6: sin(phi) * l_f = h_muscle
        sin_phi = h_muscle / fiber_length
        pennation_angle = np.arcsin(sin_phi)
        pennation_angle = np.clip(pennation_angle, PHI_MIN, PHI_MAX)

        # 3. derive l_t from g5: l_t = LMTU - cos(phi) * l_f
        tendon_length = lmtu - np.cos(pennation_angle) * fiber_length


    else: #is not stretched
        # guess by the fiber shorting
        # ----- Initial guess géométrique tenant compte de a et de l_mt -----
        # Hypothèse : déformation tendineuse croissante avec l'activation
        # (tendon ~rigide au repos, ~4% d'étirement à activation maximale)
        epsilon_f_max = 0.35
        fiber_length = l0m * (1.0 - epsilon_f_max * a) # shortening fiber

        # Petite perturbation aléatoire pour éviter de démarrer toujours
        # au même point (utile si tu fais du multi-start) 5%
        fiber_length *= (1.0 + rng.uniform(-0.005, 0.005))

        # Sécurité : on reste dans les bornes admissibles de g6
        fiber_length = np.clip(fiber_length, l_f_min, l_f_max)

        # 2. derive phi from g6: sin(phi) * l_f = h_muscle
        sin_phi = h_muscle / fiber_length
        pennation_angle = np.arcsin(sin_phi)
        pennation_angle = np.clip(pennation_angle, PHI_MIN, PHI_MAX)

        # 3. derive l_t from g5: l_t = LMTU - cos(phi) * l_f
        tendon_length = lmtu - np.cos(pennation_angle) * fiber_length

    return np.array([fiber_length, pennation_angle, tendon_length])

def get_named_params(parameters, param_labels, required):
    """
    Extract specific parameters by name from a labelled vector.

    Parameters
    ----------
    parameters : array-like
    param_labels : list of str
    required : list of str
        Names of the parameters to extract, in desired order.

    Returns
    -------
    tuple of float
        Values in the same order as `required`.

    Raises
    ------
    KeyError if any required label is missing from param_labels.
    """
    params = dict(zip(param_labels, parameters))
    missing = [k for k in required if k not in params]
    if missing:
        raise KeyError(f"Missing parameters: {missing}. "
                       f"Available: {list(params.keys())}")
    return tuple(params[k] for k in required)

def hypothetical_data_generator(skeleton_num, muscle_tendon_parameters_num, casadi_function):
    """ generate hypothetical data to make the NLP
   .
    rooted = vertcat(tendon_force,muscle_force,tendon_length, fiber_length, pennation_angle)


           Returns:
        hypothetical_data (np.ndarray): Shape (15,ntrials) including
        [ankle_torque, q_knee, q_ankle
        a_tibialis, a_soleus, a_gastrocnemius,
        fiber_length_tibialis, fiber_length_soleus, fiber_length_gastrocnemius,
        pennation_angle_tibialis, pennation_angle_soleus, pennation_angle_gastrocnemius,
        tendon_length_tibialis,tendon_length_soleus,tendon_length_gastrocnemius]
   """

    header = [
        'ankle_torque',
        'q_knee', 'q_ankle',
        'a_tibialis', 'a_soleus', 'a_gastrocnemius',
        'fiber_length_tibialis', 'fiber_length_soleus', 'fiber_length_gastrocnemius',
        'pennation_angle_tibialis', 'pennation_angle_soleus', 'pennation_angle_gastrocnemius',
        'tendon_length_tibialis', 'tendon_length_soleus', 'tendon_length_gastrocnemius'
    ]

    n_trials_fail = 0  # compteur of trials non - optimized(p > 10e-5)
    n_trials_succeeds = 0  # compteur of trials optimized(p > 10e-5)

    # ========= 1 Muscle-Tendon Architecture Equations   ========= #
    # ========= 1.1 Musculo skeletal configuration during trial(input)   ========= #
    q_knee, q_ankle, a_num = build_protocol(ramp_levels=np.array([0.1, 0.2, 0.3, 0.4, 0.9]),
                   passive_step=1.0,
                   to_radians=True)


    # ========= 2.1 data generator   ========= #
    n_trials = len(a_num)
    hypothetical_data = np.zeros((n_trials, 15))

    trials = 0

    for i in range(n_trials):  # for each measure
        # 2.1 Progression
        trials += 1
        percent = round(trials / n_trials, 2)
        print(f"Progressing: {percent * 100:.0f} %")

        # 2.2 Neuromusculoskeletal configuration
        q_num = [0, 0, 0, 0, q_knee[i], q_ankle[i]]
        musculoskeletal_states_num = q_num + list(skeleton_num)
        neuromusculoskeletal_state_num = np.concatenate([a_num[i], musculoskeletal_states_num])

        # 2.3 UMT length
        mtu_length = casadi_function['get_mtu_length'](musculoskeletal_states_num)

        # 2.4 root variables
        results = solve_all_muscles(a_num[i], mtu_length, muscle_tendon_parameters_num, casadi_function)

        x_opt_tibialis = results['tibialis']['x_opt']
        x_opt_soleus = results['soleus']['x_opt']
        x_opt_gastrocnemius = results['gastrocnemius']['x_opt']

        state_tibialis = results['tibialis']['state']
        state_soleus = results['soleus']['state']
        state_gastrocnemius = results['gastrocnemius']['state']

        if state_tibialis == 'fail' or state_soleus == 'fail' or state_gastrocnemius == 'fail':
            n_trials_fail += 1
        else:
            n_trials_succeeds += 1

            rooted_variables = np.array([x_opt_tibialis[0, 0], x_opt_soleus[0, 0], x_opt_gastrocnemius[0, 0],
                                         x_opt_tibialis[1, 0], x_opt_soleus[1, 0], x_opt_gastrocnemius[1, 0],
                                         x_opt_tibialis[2, 0], x_opt_soleus[2, 0], x_opt_gastrocnemius[2, 0]]).flatten()

            all_state = np.concatenate([neuromusculoskeletal_state_num, rooted_variables])

            ankle_torque = casadi_function['get_joint_moment'](all_state, muscle_tendon_parameters_num)

            # 2.7 extract variable
            # rooted = vertcat(tendon length,fiber length,,pennationAngle)
            tibialis_fiber_length = float(x_opt_tibialis[0])
            tibialis_pennation_angle = float(x_opt_tibialis[1])
            tibialis_tendon_length = float(x_opt_tibialis[2])

            soleus_fiber_length = float(x_opt_soleus[0])
            soleus_pennation_angle = float(x_opt_soleus[1])
            soleus_tendon_length = float(x_opt_soleus[2])

            gastrocnemius_fiber_length = float(x_opt_gastrocnemius[0])
            gastrocnemius_pennation_angle = float(x_opt_gastrocnemius[1])
            gastrocnemius_tendon_length = float(x_opt_gastrocnemius[2])

            hypothetical_data[n_trials_succeeds - 1] = [
                float(ankle_torque),
                q_knee[i], q_ankle[i],
                float(a_num[i, 0]), float(a_num[i, 1]), float(a_num[i, 2]),
                tibialis_fiber_length, soleus_fiber_length, gastrocnemius_fiber_length,
                tibialis_pennation_angle, soleus_pennation_angle, gastrocnemius_pennation_angle,
                tibialis_tendon_length, soleus_tendon_length, gastrocnemius_tendon_length
            ]

    hypothetical_data = hypothetical_data[0:n_trials_succeeds]

    hypothetical_data = hypothetical_data.T
    print('\n================== hypotetical data generator state ==================')
    print(f"fail trials: {(n_trials_fail / n_trials) * 100:.2f} %")

    return header, hypothetical_data

def build_protocol(ramp_levels=np.array([0.1, 0.2, 0.3, 0.4, 0.9]),
                   passive_step=1.0,
                   coact_jitter=0.03,
                   to_radians=True,
                   seed=0):
    """
    Construit le protocole expérimental complet en fusionnant les 3 sous-protocoles.

    Protocole 1 : Mobilisation passive
        q_knee = 0°, q_ankle ∈ [-25°, 20°] (pas = passive_step), a_num = [0, 0, 0]

    Protocole 2 : Rampes isométriques, genou tendu (q_knee = 0°)
        q_ankle ∈ {-25, -10, 0, 10, 20}
        - Dorsiflexion  : [a, 0, 0]   (tibialis seul)
        - Flexion plantaire, série A : a_sol = ramp_levels (nominal)
                                       a_gast ≈ ramp_levels (jitter ±coact_jitter)
          -> [0, a_sol, a_gast]
        - Flexion plantaire, série B : a_gast = ramp_levels (nominal)
                                       a_sol ≈ ramp_levels (jitter ±coact_jitter)
          -> [0, a_sol, a_gast]

    Protocole 3 : Rampes en dorsiflexion, genou fléchi
        q_knee = 90°, q_ankle ∈ {0, -25}, a_num = [a, 0, 0]

    Paramètres
    ----------
    ramp_levels  : array, niveaux nominaux d'activation des rampes
    passive_step : float, pas angulaire (°) du protocole passif
    coact_jitter : float, amplitude max du décalage (±) pour le muscle co-activé
                   (ex. 0.03 -> valeurs dans [nominal - 0.03, nominal + 0.03])
    to_radians   : bool, convertit les angles en radians si True
    seed         : int, graine pour la reproductibilité du jitter

    Retourne
    --------
    q_knee  : (N,) array
    q_ankle : (N,) array
    a_num   : (N, 3) array — [tibialis, soléaire, gastrocnémien]
    """

    rng = np.random.default_rng(seed)

    def _jitter(levels):
        """Décale légèrement chaque niveau autour de sa valeur nominale."""
        noise = rng.uniform(-coact_jitter, coact_jitter, size=levels.shape)
        return np.clip(levels + noise, 0.0, 1.0)

    # ---------- Protocol 1 : passive mobilization ---------- #
    qk_P1 = np.array([0.0])
    qa_P1 = np.arange(-25, 20 + passive_step, passive_step, dtype=float)
    a_P1  = np.array([[0.01, 0.01, 0.01]])

    # ---------- Protocol 2 : isometric ramp, knee extended ---------- #
    qk_P2 = np.array([0.0])
    qa_P2 = np.array([-25.0, -10.0, 0.0, 10.0, 20.0])

    # Dorsiflexion : tibialis (echo sur le tibialis)
    a_P2_DF = np.column_stack([ramp_levels,
                               np.full_like(ramp_levels, 0.01),
                               np.full_like(ramp_levels, 0.01)])

    # Plantar flexion, série A : solear nominal, gastroc jitteré (echo sur le gastroc)
    a_P2_PF_A = np.column_stack([np.full_like(ramp_levels, 0.01),
                                 ramp_levels,
                                 _jitter(ramp_levels)])

    # Plantar flexion,, série B : gastroc nominal, solear jitteré (echo sur le solear)
    a_P2_PF_B = np.column_stack([np.full_like(ramp_levels, 0.01),
                                 _jitter(ramp_levels),
                                 ramp_levels])

    a_P2 = np.vstack([a_P2_DF, a_P2_PF_A, a_P2_PF_B])   # 5 + 5 + 5 = 15 level measured

    # ---------- Protocol 3 : isometric ramp, knee flexed ---------- # (echo sur le gastroc)

    qk_P3 = np.array([90.0])
    qa_P3 = np.array([0.0, -25.0])
    a_P3  = np.column_stack([np.full_like(ramp_levels, 0.01),
                             _jitter(ramp_levels),
                             ramp_levels])


    # ---------- Product cartésien par Protocol ---------- #
    def _cartesian(qk, qa, a):
        QK, QA, A_idx = np.meshgrid(qk, qa, np.arange(len(a)), indexing='ij')
        return QK.ravel(), QA.ravel(), a[A_idx.ravel()]

    qk1, qa1, a1 = _cartesian(qk_P1, qa_P1, a_P1)
    qk2, qa2, a2 = _cartesian(qk_P2, qa_P2, a_P2)
    qk3, qa3, a3 = _cartesian(qk_P3, qa_P3, a_P3)

    # ---------- Fusion ---------- #
    q_knee  = np.concatenate([qk1, qk2, qk3])
    q_ankle = np.concatenate([qa1, qa2, qa3])
    a_num   = np.vstack([a1, a2, a3])


    if to_radians:
        q_knee  = np.deg2rad(q_knee)
        q_ankle = np.deg2rad(q_ankle)

    return q_knee, q_ankle, a_num

def simulation_(skeleton_num, muscle_tendon_parameters_num, casadi_function, time_num, q5_num, q6_num, a_tibialis_num,
                a_soleus_num, a_gastrocnemius_num):
    """ simulation_
   .
    rooted = vertcat(tendon_force,muscle_force,tendon_length, fiber_length, pennation_angle)


           Returns:
               simd_data = (np.ndarray): Shape (16,n_trials) including
               [torque,
               fiber_force_tibialis, fiber_force_soleus, fiber_force_gastrocnemius,
               fiber_active_force_tibialis, fiber_active_force_soleus, fiber_force_active_gastrocnemius,
               fiber_passive_force_tibialis, fiber_passive_force_soleus, fiber_passive_force_gastrocnemius,
               fiber_length_tibialis, fiber_length_soleus, fiber_length_gastrocnemius,
               pennation_angle_tibialis, pennation_angle_soleus, pennation_angle_gastrocnemius]
   """

    # ========= 0 set up & header & output  ========= #
    n_trials = len(q6_num)

    header_neuro_musculo_state = [
        'frame', 'time',
        'q5', 'q6',
        'a_tibialis', 'a_soleus', 'a_gastrocnemius'
    ]
    simulated_data_neuro_musculo_state = np.zeros((n_trials, len(header_neuro_musculo_state)))

    header_mtu_architecture = [
        'frame', 'time',
        'tibialis_mtu_length', 'soleus_mtu_length', 'gastrocnemius_mtu_length',
        'tibialis_fiber_length', 'soleus_fiber_length', 'gastrocnemius_fiber_length',
        'tibialis_pennation_angle', 'soleus_pennation_angle', 'gastrocnemius_pennation_angle',
        'tibialis_tendon_length', 'soleus_tendon_length', 'gastrocnemius_tendon_length'
    ]
    simulated_data_mtu_architecture = np.zeros((n_trials, len(header_mtu_architecture)))

    header_mtu_forces = ['frame', 'time',
                         'torque_simulated',
                         'tendon_force_tibialis', 'tendon_force_soleus', 'tendon_force_gastrocnemius',
                         'fiber_force_tibialis', 'fiber_force_soleus', 'fiber_force_gastrocnemius',
                         'fiber_active_force_tibialis', 'fiber_active_force_soleus', 'fiber_active_force_gastrocnemius',
                         'fiber_passive_force_tibialis', 'fiber_passive_force_soleus',
                         'fiber_passive_force_gastrocnemius']
    simulated_data_mtu_forces = np.zeros((n_trials, len(header_mtu_forces)))

    n_trials_fail = 0  # compteur of trials non - optimized(p > 10e-5)
    n_trials_succeds = 0  # compteur of trials optimized(p > 10e-5)

    # ========= 1 Muscle-Tendon Architecture Equations   ========= #
    # ========= 1.1 Musculo skeletical configuration during trial(input)   ========= #

    q5_num = np.deg2rad(q5_num)
    q6_num = np.deg2rad(q6_num)

    # ========= 1.2 Neuronal activation(input)   ========= #
    # Muscle activation [a_tibialis_num, a_soleus_num, a_gastrocnemius_num]

    # ========= 1.3 Muscle Tendon Parameters(ℓom, φo, Fom, ℓst)(input)   ========= #
    muscle_tendon_parameters_ta = np.array(muscle_tendon_parameters_num)[[0, 3, 6, 9]]
    muscle_tendon_parameters_sol = np.array(muscle_tendon_parameters_num)[[1, 4, 7, 10]]
    muscle_tendon_parameters_gast = np.array(muscle_tendon_parameters_num)[[2, 5, 8, 11]]

    # ========= 2.1 data generator   ========= #

    trials = 0

    for i in range(n_trials):  # for each frame
        # 2.1 Progression
        trials += 1
        percent = round(trials / n_trials, 2)
        print(f"Progressing: {percent * 100:.0f} %")

        # 2.2 Neuromusculoskeletal configuration
        a_num = [a_tibialis_num[i], a_soleus_num[i], a_gastrocnemius_num[i]]
        q_num = [0, 0, 0, 0, q5_num[i], q6_num[i]]

        musculoskeletal_states_num = q_num + list(skeleton_num)
        neuromusculoskeletal_state_num = np.concatenate([a_num, musculoskeletal_states_num])

        # 2.3 UMT length
        mtu_length = casadi_function['get_mtu_length'](musculoskeletal_states_num)
        tibialis_mtu_length = float(mtu_length[0])
        soleus_mtu_length = float(mtu_length[1])
        gastrocnemius_mtu_length = float(mtu_length[2])

        # 2.4 muscle activities
        a_tibialis = a_num[0]
        a_soleus = a_num[1]
        a_gastrocnemius = a_num[2]

        # 2.5 Solver
        x_opt_tibialis, state_tibialis = root_muscle_dynamics(a_tibialis,
                                                              tibialis_mtu_length,
                                                              muscle_tendon_parameters_ta,
                                                              'tibialis',
                                                              casadi_function)

        x_opt_soleus, state_soleus = root_muscle_dynamics(a_soleus,
                                                          soleus_mtu_length,
                                                          muscle_tendon_parameters_sol,
                                                          'soleus',
                                                          casadi_function)

        x_opt_gastrocnemius, state_gastrocnemius = root_muscle_dynamics(a_gastrocnemius,
                                                                        gastrocnemius_mtu_length,
                                                                        muscle_tendon_parameters_gast,
                                                                        'gastrocnemius',
                                                                        casadi_function)

        if state_tibialis == 'fail' or state_soleus == 'fail' or state_tibialis == 'fail' or state_gastrocnemius == 'fail':
            n_trials_fail += 1
        else:
            n_trials_succeds += 1

            rooted_fiber_length = np.array(
                [x_opt_tibialis[0, 0], x_opt_soleus[0, 0], x_opt_gastrocnemius[0, 0]]).flatten()
            rooted_pennation_angle = np.array(
                [x_opt_tibialis[1, 0], x_opt_soleus[1, 0], x_opt_gastrocnemius[1, 0]]).flatten()
            rooted_tendon_length = np.array(
                [x_opt_tibialis[2, 0], x_opt_soleus[2, 0], x_opt_gastrocnemius[2, 0]]).flatten()

            rooted_variables = np.concatenate([rooted_fiber_length, rooted_pennation_angle, rooted_tendon_length])

            all_state = np.concatenate([neuromusculoskeletal_state_num, rooted_variables])

            # 2.7 extract variable
            # about architecture
            # rooted = vertcat(tendon length,fiber lenght,,pennationAngle)
            tibialis_fiber_length = float(x_opt_tibialis[0])
            tibialis_pennation_angle = x_opt_tibialis[1]
            tibialis_pennation_angle = float(np.rad2deg(tibialis_pennation_angle))
            tibialis_tendon_length = float(x_opt_tibialis[2])

            soleus_fiber_length = float(x_opt_soleus[0])
            soleus_pennation_angle = x_opt_soleus[1]
            soleus_pennation_angle = float(np.rad2deg(soleus_pennation_angle))
            soleus_tendon_length = float(x_opt_soleus[2])

            gastrocnemius_fiber_length = float(x_opt_gastrocnemius[0])
            gastrocnemius_pennation_angle = x_opt_gastrocnemius[1]
            gastrocnemius_pennation_angle = float(np.rad2deg(gastrocnemius_pennation_angle))
            gastrocnemius_tendon_length = float(x_opt_gastrocnemius[2])

            # about forces
            ankle_torque = casadi_function['get_joint_moment'](all_state, muscle_tendon_parameters_num)
            torque_simulated = float(ankle_torque)

            tendon_force = casadi_function['get_tendon_force'](all_state, muscle_tendon_parameters_num)
            tendon_force_tibialis = float(tendon_force[0])
            tendon_force_soleus = float(tendon_force[1])
            tendon_force_gastrocnemius = float(tendon_force[2])

            fiber_force = casadi_function['get_fiber_force'](all_state, muscle_tendon_parameters_num)
            fiber_force_tibialis = float(fiber_force[0])
            fiber_force_soleus = float(fiber_force[1])
            fiber_force_gastrocnemius = float(fiber_force[2])

            fiber_active_force = casadi_function['get_fiber_active_force'](all_state, muscle_tendon_parameters_num)
            fiber_active_force_tibialis = float(fiber_active_force[0])
            fiber_active_force_soleus = float(fiber_active_force[1])
            fiber_active_force_gastrocnemius = float(fiber_active_force[2])

            fiber_passive_force = casadi_function['get_fiber_passive_force'](all_state, muscle_tendon_parameters_num)
            fiber_passive_force_tibialis = float(fiber_passive_force[0])
            fiber_passive_force_soleus = float(fiber_passive_force[1])
            fiber_passive_force_gastrocnemius = float(fiber_passive_force[2])

            # 2.output
            # header_neuro_musculo_state = [
            #     'frame', 'time',
            #     'q5', 'q6',
            #     'a_tibialis', 'a_soleus', 'a_gastrocnemius'
            # ]

            simulated_data_neuro_musculo_state[i, :] = [
                trials, time_num[i],
                np.rad2deg(q_num[4]), np.rad2deg(q_num[5]),
                a_num[0] * 100, a_num[1] * 100, a_num[2] * 100
            ]
            # header_mtu_architecture = [
            #     'frame', 'time',
            #     'tibialis_mtu_length', 'soleus_mtu_length', 'gastrocnemius_mtu_length',
            #     'tibialis_fiber_length', 'soleus_fiber_length', 'gastrocnemius_fiber_length',
            #     'tibialis_pennation_angle', 'soleus_pennation_angle','gastrocnemius_pennation_angle',
            #     'tibialis_tendon_length', 'soleus_tendon_length', 'gastrocnemius_tendon_length'
            # ]
            simulated_data_mtu_architecture[i, :] = [
                trials, time_num[i],
                tibialis_mtu_length, soleus_mtu_length, gastrocnemius_mtu_length,
                tibialis_fiber_length, soleus_fiber_length, gastrocnemius_fiber_length,
                tibialis_pennation_angle, soleus_pennation_angle, gastrocnemius_pennation_angle,
                tibialis_tendon_length, soleus_tendon_length, gastrocnemius_tendon_length
            ]

            # header_mtu_forces = [
            #     'frame', 'time',
            #     'torque_simulated',
            #     'tendon_force_tibialis', 'tendon_force_soleus', 'tendon_force_gastrocnemius',
            #     'fiber_force_tibialis', 'fiber_force_soleus', 'fiber_force_gastrocnemius',
            #     'fiber_active_force_tibialis', 'fiber_active_force_soleus', 'fiber_active_force_gastrocnemius',
            #     'fiber_passive_force_tibialis', 'fiber_passive_force_soleus', 'fiber_passive_force_gastrocnemius'
            # ]
            simulated_data_mtu_forces[i, :] = [
                trials, time_num[i],
                torque_simulated,
                tendon_force_tibialis, tendon_force_soleus, tendon_force_gastrocnemius,
                fiber_force_tibialis, fiber_force_soleus, fiber_force_gastrocnemius,
                fiber_active_force_tibialis, fiber_active_force_soleus, fiber_active_force_gastrocnemius,
                fiber_passive_force_tibialis, fiber_passive_force_soleus, fiber_passive_force_gastrocnemius
            ]

    print('\n================== simulation exit ==================')
    print(f"fail trials: {(n_trials_fail / n_trials) * 100:.2f} %")

    return header_neuro_musculo_state, simulated_data_neuro_musculo_state, header_mtu_architecture, simulated_data_mtu_architecture, header_mtu_forces, simulated_data_mtu_forces

def add_noise(data, sigma_torque=0.1, sigma_length=0.002, sigma_angle_deg=3.0,
              add_internal_noise=False, rng=None):
    """
    Add Gaussian sensor noise to hypothetical_data IN PLACE.

    Noise is applied to measurable quantities only (torque, joint angles).
    Internal Hill-model variables (fiber length, pennation, tendon length,
    activations) are left untouched by default — set add_internal_noise=True
    to also perturb them, e.g. to simulate ultrasound measurement noise on
    fiber length and pennation.

    Parameters
    ----------
    data : np.ndarray, shape (15, n_trials)
        Hypothetical data matrix. See module docstring for row layout.
        Modified in place.
    sigma_torque : float, default 0.1
        Std of torque noise [N.m].
    sigma_length : float, default 0.002
        Std of length noise [m] (= 0.2 cm).
    sigma_angle_deg : float, default 3.0
        Std of angle noise [degrees]. Internally converted to radians.
    add_internal_noise : bool, default False
        If True, also adds noise to fiber length, pennation angle, and
        tendon length (typical for ultrasound-based estimates).
    rng : np.random.Generator, optional
        For reproducibility. If None, uses np.random.default_rng().

    Returns
    -------
    np.ndarray
        The same data array (modified in place), returned for chaining.
    """
    # Indices des lignes dans la matrice (15, n_trials)
    ROW_TORQUE = 0  # 1 ligne   : couple [N.m]
    ROW_Q = slice(1, 3)  # 2 lignes  : angles articulaires [rad]
    ROW_ACTIVATION = slice(3, 6)  # 3 lignes  : activations [0, 1]
    ROW_FIBER_LENGTH = slice(6, 9)  # 3 lignes  : longueurs fibres [m]
    ROW_PENNATION = slice(9, 12)  # 3 lignes : angles pennation [rad]
    ROW_TENDON_LENGTH = slice(12, 15)  # 3 lignes : longueurs tendon [m]


    rng = rng if rng is not None else np.random.default_rng()
    sigma_angle_rad = np.deg2rad(sigma_angle_deg)
    n_trials = data.shape[1]

    # === Bruit sur les mesures (capteur) === #
    # Couple : 1 ligne
    data[ROW_TORQUE, :] += rng.normal(0.0, sigma_torque, size=n_trials)

    # Angles articulaires : 2 lignes
    data[ROW_Q, :] += rng.normal(0.0, sigma_angle_rad/8, size=(2, n_trials))

    # === Bruit sur les estimations échographiques (optionnel) === #
    if add_internal_noise:
        # Longueurs de fibres : 3 lignes
        data[ROW_FIBER_LENGTH, :] += rng.normal(0.0, sigma_length, size=(3, n_trials))

        # Angles de pennation : 3 lignes
        data[ROW_PENNATION, :] += rng.normal(0.0, sigma_angle_rad, size=(3, n_trials))

        # Longueurs de tendon : 3 lignes
        data[ROW_TENDON_LENGTH, :] += rng.normal(0.0, sigma_length, size=(3, n_trials))

    return data


def nlp_verification(data, initial_guess, lower_band, upper_band, skeleton_num,
                     muscle_tendon_parameters_num, unknown_parameters, casadi_function):
    # 1. Construire w0 aux vrais paramètres ET vrais états
    #    (les états qui ont servi à générer les données)
    w0_true = np.concatenate([
        muscle_tendon_parameters_num,  # les 12 vrais params
        # puis pour chaque trial, les états utilisés pour générer les mesures :
        # soit mesured_fiber_length, mesured_pennation_angle, mesured_tendon_length
        # concaténés dans le bon ordre (fl, pa, tl) pour chaque trial
    ])

    """
    Identify muscle-tendon parameters (ℓom, φo, Fom, ℓst) via an NLP.

    Args:
        data (np.ndarray): Shape (15, n_trials).
            Rows:
              [0]     : measured joint torque (N.m)
              [1:3]   : joint angles q (rad)
              [3:6]   : muscle activations [0, 1]
              [6:9]   : measured fiber lengths (m)
              [9:12]  : measured pennation angles (rad)
              [12:15] : measured tendon lengths (m)
        initial_guess, lower_band, upper_band : muscle-tendon parameters (size 12).
        skeleton_num : skeleton geometry (.osim).
        muscle_tendon_parameters_num : reference values used for comparison.
        unknown_parameters : SX vector of the 12 unknowns.
        casadi_function : dict of CasADi functions.

    Returns:
        xopt (np.ndarray): Shape (12,) — estimated (ℓom, φo, Fom, ℓst).
    """

    # ============ Sanity checks ============ #
    assert data.shape[0] == 15, f"data must be (15, n_trials), got {data.shape}"

    lower_band = np.asarray(lower_band, dtype=float)
    upper_band = np.asarray(upper_band, dtype=float)
    initial_guess = np.asarray(initial_guess, dtype=float)

    assert np.all(lower_band <= upper_band), \
        f"lower_band > upper_band at: {np.where(lower_band > upper_band)[0]}"

    n_trials = data.shape[1]
    n_muscle = 3

    # ============ NLP set up ============ #
    # Decision variables, initial guess, and bounds
    w, w0, lbw, ubw = [], [], [], []
    j = 0
    # Constraints and their bounds
    g, lbg, ubg = [], [], []

    up = unknown_parameters

    # Per-trial residuals (kept for optional post-processing)
    e_torque, e_fiber, e_pennation = [], [], []

    # Unknown muscle-tendon parameters (12 values)
    w += [up]
    w0 += list(initial_guess)
    lbw += list(lower_band)
    ubw += list(upper_band)

    # Weights in the cost function
    w_torque = 1  # N.m
    w_length = 0.005  # mm
    w_angle = (1 / 180) * np.pi  # rad

    for trial in range(n_trials):
        # --- Extract data for this trial (one column = one trial) ---
        data_trial = data[:, trial]

        measured_torque = data_trial[0]  # N.m
        q_trial = [0, 0, 0, 0] + list(data_trial[1:3])  # rad
        a_trial = data_trial[3:6]  # [0, 1]
        measured_fiber_length = data_trial[6:9]  # m
        measured_pennation_angle = data_trial[9:12]  # rad
        measured_tendon_length = data_trial[12:15]  # m

        musculoskeletal_states_trial = q_trial + list(skeleton_num)
        neuromusculoskeletal_state_trial = np.concatenate(
            [a_trial, musculoskeletal_states_trial]
        )

        # MTU (muscle-tendon unit) length
        mtu_length = casadi_function['get_mtu_length'](musculoskeletal_states_trial)

        # --- Trial-specific decision variables ---
        tendon_length_k = SX.sym(f"Tendon_Length_{trial + 1}", n_muscle)
        fiber_length_k = SX.sym(f"Fiber_Length_{trial + 1}", n_muscle)
        pennation_angle_k = SX.sym(f"Pennation_Angle_{trial + 1}", n_muscle)

        # --- Physiological bounds (robust to sign) ---
        fl_meas = np.abs(measured_fiber_length)
        tl_meas = np.abs(measured_tendon_length)

        # Fiber length: positive, around the measured value
        lb_fl = fl_meas - fl_meas * 0.3
        ub_fl = fl_meas + fl_meas * 0.3

        # Pennation: physiological bounds, sign is free
        lb_pa = np.full(n_muscle, -np.pi / 2 + 1e-3)
        ub_pa = np.full(n_muscle, np.pi / 2 - 1e-3)

        # Tendon length: positive, tight range
        lb_tl = tl_meas * 0.05
        ub_tl = tl_meas * 1.05

        lb_block = np.concatenate([lb_fl, lb_pa, lb_tl])
        ub_block = np.concatenate([ub_fl, ub_pa, ub_tl])

        # Initial guess, with lengths forced positive
        w0_k = np.concatenate([fl_meas, measured_pennation_angle, tl_meas])
        w_k = vertcat(fiber_length_k, pennation_angle_k, tendon_length_k)

        # Safety: check bounds consistency and that w0 lies within them
        assert np.all(lb_block <= ub_block), \
            f"Inverted bounds at trial {trial}"
        assert np.all((w0_k >= lb_block - 1e-9) & (w0_k <= ub_block + 1e-9)), \
            f"w0 out of bounds at trial {trial}: " \
            f"w0_k={w0_k}, lb={lb_block}, ub={ub_block}"

        w += [w_k]
        w0 += list(w0_k)
        lbw += list(lb_block)
        ubw += list(ub_block)

        # --- Equilibrium constraints (9 per trial) ---
        k = vertcat(a_trial, mtu_length, up)
        constraints = casadi_function['equilibrium_error_all_muscle'](w_k, k)
        g += [constraints]
        lbg += [0] * 9
        ubg += [0] * 9

        # --- Torque simulation ---
        all_states = vertcat(SX(neuromusculoskeletal_state_trial.tolist()), w_k)
        torque_simulated = casadi_function['get_joint_moment'](all_states, unknown_parameters)

        # --- Residuals ---
        e_torque_trials = measured_torque - torque_simulated
        e_fiber_trials = measured_fiber_length - fiber_length_k
        e_pennation_trials = measured_pennation_angle - pennation_angle_k

        j += w_torque * e_torque_trials ** 2
        j += sum1(w_length * e_fiber_trials ** 2)
        j += sum1(w_angle * e_pennation_trials ** 2)

        e_torque.append(e_torque_trials)
        e_fiber.append(e_fiber_trials)
        e_pennation.append(e_pennation_trials)

    # ============ Assembly ============ #
    w = vertcat(*w)
    g = vertcat(*g)
    print("J shape:", j.shape)
    print(f"n_trials : {n_trials} | size(w) : {w.shape} | size(g) : {g.shape}")

    w0 = np.array(w0, dtype=float)
    lbw = np.array(lbw, dtype=float)
    ubw = np.array(ubw, dtype=float)

    # Global checks
    if np.any(np.isnan(w0)):
        print('NaNs in w0 at:', np.where(np.isnan(w0)))
        return None
    if np.any(np.isinf(w0)):
        print('Infs in w0 at:', np.where(np.isinf(w0)))
        return None
    bad = np.where(lbw > ubw)[0]
    if len(bad):
        print(f'lbw > ubw at indices: {bad[:20]}...')
        return None

    print('w0 is valid')

    # ============ NLP solver ============ #
    nlp = {'x': w, 'f': j, 'g': g}
    solver = nlpsol('solver', 'ipopt', nlp)
    print(solver)
    w0_true_list = [np.asarray(muscle_tendon_parameters_num, dtype=float).flatten()]

    # 3. Évaluer
    for trial in range(n_trials):
        data_trial = data[:, trial]
        measured_fiber_length = data_trial[6:9]
        measured_pennation_angle = data_trial[9:12]
        measured_tendon_length = data_trial[12:15]

        # ATTENTION : même ordre que w_k = vertcat(fiber, pennation, tendon)
        # et même traitement (valeur absolue sur les longueurs comme dans w0_k)
        fl = np.abs(measured_fiber_length)
        pa = np.asarray(measured_pennation_angle, dtype=float).flatten()
        tl = np.abs(measured_tendon_length)

        w0_true_list.append(np.concatenate([fl, pa, tl]))

    w0_true = np.concatenate(w0_true_list)

    # Vérification de la taille
    expected_size = 12 + 9 * n_trials
    assert w0_true.shape[0] == expected_size, \
        f"Expected {expected_size}, got {w0_true.shape[0]}"
    print(f"w0_true size: {w0_true.shape[0]} (expected {expected_size})")

    # --- Évaluer f et g à ce point ---
    f_eval = Function('f_eval', [w], [j])
    g_eval = Function('g_eval', [w], [g])

    f_val = float(f_eval(w0_true))
    g_val = np.array(g_eval(w0_true)).flatten()

    print(f"\nCoût aux vrais paramètres/états : {f_val:.6e}")
    print(f"Max |g| aux vrais params/états  : {np.max(np.abs(g_val)):.6e}")
    print(f"Nb contraintes violées > 1e-6   : {np.sum(np.abs(g_val) > 1e-6)}")
    print(f"Nb contraintes violées > 1e-3   : {np.sum(np.abs(g_val) > 1e-3)}")

    # Setup du trial 0
    trial = 0
    data_trial = data[:, trial]
    a_trial = data_trial[3:6]
    q_trial = [0, 0, 0, 0] + list(data_trial[1:3])
    musculoskeletal_states_trial = q_trial + list(skeleton_num)
    mtu_length = casadi_function['get_mtu_length'](musculoskeletal_states_trial)

    fl_meas = np.abs(data_trial[6:9])
    pa_meas = data_trial[9:12]
    tl_meas = np.abs(data_trial[12:15])
    lom = muscle_tendon_parameters_num[0:3]
    lst = muscle_tendon_parameters_num[9:12]

    k = vertcat(a_trial, mtu_length, muscle_tendon_parameters_num)
    f_eq = casadi_function['equilibrium_error_all_muscle']

    # --- TEST 1 : comme dans ton code actuel ---
    w_k_test1 = np.concatenate([fl_meas, pa_meas, tl_meas])
    g1 = np.array(f_eq(w_k_test1, k)).flatten()
    print(f"Test 1 (ordre fl,pa,tl, m)           : max|g| = {np.max(np.abs(g1)):.3e}")

    # --- TEST 2 : fiber et tendon normalisées ---
    w_k_test2 = np.concatenate([fl_meas / lom, pa_meas, tl_meas / lst])
    g2 = np.array(f_eq(w_k_test2, k)).flatten()
    print(f"Test 2 (fl/lom, pa, tl/lst)          : max|g| = {np.max(np.abs(g2)):.3e}")

    # --- TEST 3 : que la fiber normalisée ---
    w_k_test3 = np.concatenate([fl_meas / lom, pa_meas, tl_meas])
    g3 = np.array(f_eq(w_k_test3, k)).flatten()
    print(f"Test 3 (fl/lom, pa, tl en m)         : max|g| = {np.max(np.abs(g3)):.3e}")

    # --- TEST 4 : ordre permuté fiber, tendon, pennation ---
    w_k_test4 = np.concatenate([fl_meas, tl_meas, pa_meas])
    g4 = np.array(f_eq(w_k_test4, k)).flatten()
    print(f"Test 4 (ordre fl,tl,pa en m)         : max|g| = {np.max(np.abs(g4)):.3e}")

    # --- TEST 5 : ordre tendon, fiber, pennation ---
    w_k_test5 = np.concatenate([tl_meas, fl_meas, pa_meas])
    g5 = np.array(f_eq(w_k_test5, k)).flatten()
    print(f"Test 5 (ordre tl,fl,pa en m)         : max|g| = {np.max(np.abs(g5)):.3e}")

    # --- TEST 6 : pennation avec variable normalisée cos(phi) ---
    # (au cas où ce serait déjà cos(phi) et pas phi)
    w_k_test6 = np.concatenate([fl_meas, np.cos(pa_meas), tl_meas])
    g6 = np.array(f_eq(w_k_test6, k)).flatten()
    print(f"Test 6 (fl, cos(pa), tl)             : max|g| = {np.max(np.abs(g6)):.3e}")

    f_eq = casadi_function['equilibrium_error_all_muscle']
    print(f"Entrée 0 '{f_eq.name_in(0)}' : shape {f_eq.size_in(0)}")
    print(f"Entrée 1 '{f_eq.name_in(1)}' : shape {f_eq.size_in(1)}")

    # Setup trial 0
    trial = 0
    data_trial = data[:, trial]
    a_trial = data_trial[3:6]
    q_trial = [0, 0, 0, 0] + list(data_trial[1:3])
    musculoskeletal_states_trial = q_trial + list(skeleton_num)
    mtu_length = casadi_function['get_mtu_length'](musculoskeletal_states_trial)

    # IMPORTANT : reconvertir en objets numériques concrets
    mtu_length_num = np.array(mtu_length).flatten()
    p_num = np.concatenate([
        np.asarray(a_trial, dtype=float).flatten(),
        mtu_length_num,
        np.asarray(muscle_tendon_parameters_num, dtype=float).flatten()
    ])
    print(f"p_num shape : {p_num.shape} (doit être 18)")

    # Construire le rootfinder avec la bonne signature
    x_sym = SX.sym('x', 9)
    p_sym = SX.sym('p', 18)
    residual = casadi_function['equilibrium_error_all_muscle'](x_sym, p_sym)

    # Attention : 'g' doit dépendre de x ET p, avec p comme paramètre du rootfinder
    rf_problem = {'x': x_sym, 'p': p_sym, 'g': residual}
    rf = rootfinder('rf', 'newton', rf_problem)

    # Initial guess = mesures
    w0_guess = np.concatenate([
        np.abs(data_trial[6:9]),
        data_trial[9:12],
        np.abs(data_trial[12:15])
    ])

    # Pour un modèle de Hill standard : mtu = tendon + fiber * cos(pennation)
    trial = 0
    data_trial = data[:, trial]

    fiber_stored = np.abs(data_trial[6:9])
    pennation_stored = data_trial[9:12]
    tendon_stored = np.abs(data_trial[12:15])

    mtu_expected = tendon_stored + fiber_stored * np.cos(pennation_stored)

    # mtu_length actuel
    q_trial = [0, 0, 0, 0] + list(data_trial[1:3])
    musculoskeletal_states_trial = q_trial + list(skeleton_num)
    mtu_current = np.array(
        casadi_function['get_mtu_length'](musculoskeletal_states_trial)
    ).flatten()

    print(f"mtu attendue (depuis états stockés) : {mtu_expected}")
    print(f"mtu actuelle (get_mtu_length)       : {mtu_current}")
    print(f"Écart                               : {mtu_current - mtu_expected}")

    # Appel : rootfinder prend (x0, p)
    try:
        w_k_fresh = np.array(rf(w0_guess, p_num)).flatten()
        print("\nRootfinder OK")
        print(f"États régénérés (trial {trial}):")
        print(f"  fiber     : {w_k_fresh[0:3]}")
        print(f"  pennation : {w_k_fresh[3:6]}")
        print(f"  tendon    : {w_k_fresh[6:9]}")

        print(f"\nÉtats stockés dans data :")
        print(f"  fiber     : {np.abs(data_trial[6:9])}")
        print(f"  pennation : {data_trial[9:12]}")
        print(f"  tendon    : {np.abs(data_trial[12:15])}")

        print(f"\nÉcart absolu :")
        print(f"  fiber     : {w_k_fresh[0:3] - np.abs(data_trial[6:9])}")
        print(f"  pennation : {w_k_fresh[3:6] - data_trial[9:12]}")
        print(f"  tendon    : {w_k_fresh[6:9] - np.abs(data_trial[12:15])}")

        # Vérifier que g est bien ≈ 0 sur les états régénérés
        g_fresh = np.array(
            casadi_function['equilibrium_error_all_muscle'](w_k_fresh, p_num)
        ).flatten()
        print(f"\n|g| max aux états régénérés : {np.max(np.abs(g_fresh)):.3e}")

        # Et à nouveau avec les états stockés
        w_stored = np.concatenate([
            np.abs(data_trial[6:9]),
            data_trial[9:12],
            np.abs(data_trial[12:15])
        ])
        g_stored = np.array(
            casadi_function['equilibrium_error_all_muscle'](w_stored, p_num)
        ).flatten()
        print(f"|g| max aux états stockés   : {np.max(np.abs(g_stored)):.3e}")

    except Exception as e:
        print(f"Rootfinder a échoué : {e}")


def optimization_nlp(data, initial_guess, lower_band, upper_band, skeleton_num,
                     muscle_tendon_parameters_num, unknown_parameters,
                     casadi_function, param_index):
    """
    Identify the muscle-tendon parameters declared 'sym' in param_config,
    via an NLP. The set and order of optimized parameters are driven by
    param_index (returned by get_initial_guess).

    Args:
        data (np.ndarray): Shape (15, n_trials).
            Rows:
              [0]     : measured joint torque (N.m)
              [1:3]   : joint angles q (rad)
              [3:6]   : muscle activations [0, 1]
              [6:9]   : measured fiber lengths (m)
              [9:12]  : measured pennation angles (rad)
              [12:15] : measured tendon lengths (m)
        initial_guess, lower_band, upper_band : optimized MTU parameters,
            size n_opt = 3 * n_sym, ordered as param_index.
        skeleton_num : skeleton geometry (.osim).
        muscle_tendon_parameters_num : reference values (size n_opt), same
            order as the optimized vector (only the 'sym' parameters).
        unknown_parameters : SX vector of the n_opt unknowns (only 'sym').
        casadi_function : dict of CasADi functions (fixed params already
            hard-coded inside the generated equations).
        param_index (dict): {param_name: {muscle: idx}} mapping each optimized
            parameter/muscle to its index in the vector. Drives reporting.

    Returns:
        param_opt (np.ndarray): Shape (n_opt,) — estimated parameters,
            same order as initial_guess / param_index.
    """

    # ============ Sanity checks ============ #
    assert data.shape[0] == 15, f"data must be (15, n_trials), got {data.shape}"

    initial_guess = np.asarray(initial_guess, dtype=float)
    lower_band = np.asarray(lower_band, dtype=float)
    upper_band = np.asarray(upper_band, dtype=float)

    # Taille réelle du vecteur optimisé (ne dépend plus du nombre 12 codé en dur)
    n_opt = unknown_parameters.shape[0]
    assert initial_guess.shape[0] == n_opt, \
        f"initial_guess ({initial_guess.shape[0]}) != unknown_parameters ({n_opt})"
    assert lower_band.shape[0] == n_opt and upper_band.shape[0] == n_opt, \
        "lower_band/upper_band incohérents avec unknown_parameters"

    assert np.all(lower_band <= upper_band), \
        f"lower_band > upper_band at: {np.where(lower_band > upper_band)[0]}"

    n_trials = data.shape[1]
    n_muscle = 3

    # ============ [FIX 3] Scaling des paramètres MTU ============ #
    # On optimise des ratios O(1) au lieu des valeurs physiques brutes.
    scale = initial_guess.copy()
    assert np.all(scale > 0), \
        "initial_guess doit être strictement positif pour servir d'échelle"

    up = unknown_parameters            # variable IPOPT (sans dimension, ~ O(1))
    up_phys = up * scale               # valeurs physiques (uniquement les 'sym')

    # ============ NLP set up ============ #
    w, w0, lbw, ubw = [], [], [], []
    g, lbg, ubg = [], [], []
    j = SX(0)
    j_torque = SX(0)
    j_fascicle_length = SX(0)
    j_pennation_angle = SX(0)

    e_torque, e_fiber, e_pennation = [], [], []

    # Unknown MTU parameters (n_opt values), EN UNITÉS SCALÉES
    w += [up]
    w0 += [1.0] * n_opt                       # initial_guess / scale = 1
    lbw += list(lower_band / scale)
    ubw += list(upper_band / scale)

    # ============ [FIX 4] Poids du coût = 1 / sigma**2 ============ #
    sigma_torque = 50                                 # N.m
    sigma_fiber = 0.01                                # ~5 mm en m
    sigma_penn = np.deg2rad(5)                         # ~5° en rad

    for trial in range(n_trials):
        data_trial = data[:, trial]

        measured_torque = data_trial[0]                 # N.m
        q_trial = [0, 0, 0, 0] + list(data_trial[1:3])  # rad
        a_trial = data_trial[3:6]                        # [0, 1]
        measured_fiber_length = data_trial[6:9]          # m
        measured_pennation_angle = data_trial[9:12]      # rad
        measured_tendon_length = data_trial[12:15]       # m

        musculoskeletal_states_trial = q_trial + list(skeleton_num)
        neuromusculoskeletal_state_trial = np.concatenate(
            [a_trial, musculoskeletal_states_trial]
        )

        mtu_length = casadi_function['get_mtu_length'](musculoskeletal_states_trial)

        tendon_length_k = SX.sym(f"Tendon_Length_{trial}", n_muscle)
        fiber_length_k = SX.sym(f"Fiber_Length_{trial}", n_muscle)
        pennation_angle_k = SX.sym(f"Pennation_Angle_{trial}", n_muscle)

        fl_meas = np.abs(measured_fiber_length)
        tl_meas = np.abs(measured_tendon_length)
        pa_meas = np.abs(measured_pennation_angle)

        lb_fl = fl_meas - fl_meas * 0.05
        ub_fl = fl_meas + fl_meas * 0.05

        lb_pa = pa_meas - np.deg2rad(2)
        ub_pa = pa_meas + np.deg2rad(2)
        lb_pa = np.clip(lb_pa, 0 + epsilon, np.pi - epsilon)
        ub_pa = np.clip(ub_pa, 0 + epsilon, np.pi - epsilon)

        lb_tl = tl_meas - tl_meas * 0.005
        ub_tl = tl_meas + tl_meas * 0.005

        lb_block = np.concatenate([lb_fl, lb_pa, lb_tl])
        ub_block = np.concatenate([ub_fl, ub_pa, ub_tl])

        w0_k = np.concatenate([fl_meas, pa_meas, tl_meas])
        w_k = vertcat(fiber_length_k, pennation_angle_k, tendon_length_k)

        assert np.all(lb_block <= ub_block), \
            f"Inverted bounds at trial {trial}"
        assert np.all((w0_k >= lb_block - 1e-9) & (w0_k <= ub_block + 1e-9)), \
            f"w0 out of bounds at trial {trial}: " \
            f"w0_k={w0_k}, lb={lb_block}, ub={ub_block}"

        w += [w_k]
        w0 += list(w0_k)
        lbw += list(lb_block)
        ubw += list(ub_block)

        # --- Equilibrium constraints (9 per trial) --- #
        # up_phys ne contient que les 'sym' ; les 'fixed' sont déjà codés en
        # dur dans les équations générées par get_model_equation.
        k = vertcat(a_trial, mtu_length, up_phys)
        constraints = casadi_function['equilibrium_error_all_muscle'](w_k, k)
        g += [constraints]
        lbg += [0] * 9
        ubg += [0] * 9

        # --- Torque simulation --- #
        all_states = vertcat(SX(neuromusculoskeletal_state_trial.tolist()), w_k)
        torque_simulated = casadi_function['get_joint_moment'](all_states, up_phys)

        # --- Residuals --- #
        e_torque_trials = measured_torque - torque_simulated
        e_fiber_trials = measured_fiber_length - fiber_length_k
        e_pennation_trials = measured_pennation_angle - pennation_angle_k

        j_torque += (e_torque_trials ** 2) / sigma_torque
        j_fascicle_length += sum1((e_fiber_trials ** 2) / sigma_fiber )
        """
        j_pennation_angle += sum1((e_pennation_trials ** 2) / sigma_penn ** 2)
        """

        e_torque.append(e_torque_trials)
        e_fiber.append(e_fiber_trials)
        e_pennation.append(e_pennation_trials)

    j = j_torque + j_fascicle_length + j_pennation_angle
    n_element = 4 # number of element compared by trials (1 torque, 3 fibers= 4) or (1 torque, 3 fibers + 3 pennation angle = 7)
    n_comparison = n_element * n_trials

    j =  j / n_comparison

    # ============ Assembly ============ #
    w = vertcat(*w)
    g = vertcat(*g)
    print("J shape:", j.shape)
    print(f"n_trials : {n_trials} | size(w) : {w.shape} | size(g) : {g.shape}")

    w0 = np.array(w0, dtype=float)
    lbw = np.array(lbw, dtype=float)
    ubw = np.array(ubw, dtype=float)
    lbg = np.array(lbg, dtype=float)
    ubg = np.array(ubg, dtype=float)

    if np.any(np.isnan(w0)):
        print('NaNs in w0 at:', np.where(np.isnan(w0)))
        return None
    if np.any(np.isinf(w0)):
        print('Infs in w0 at:', np.where(np.isinf(w0)))
        return None
    bad = np.where(lbw > ubw)[0]
    if len(bad):
        print(f'lbw > ubw at indices: {bad[:20]}...')
        return None

    print('w0 is valid')

    # ============ NLP solver ============ #

    from mtu_live_callback import run_with_live_plot

    grad_f_func = Function("grad_f", [w], [gradient(j, w)])
    grad_g_func = Function("grad_g", [w], [jacobian(g, w).T])

    opts_ipopt = {
        "ipopt.max_iter": 1000,
        "ipopt.tol": 1e-4,
        "ipopt.print_info_string": "yes",
        "ipopt.linear_solver": "mumps",
        "ipopt.hessian_approximation": "exact",
    }

    sol = run_with_live_plot(
        nlp={'x': w, 'f': j, 'g': g},
        x0=w0, lbx=lbw, ubx=ubw, lbg=lbg, ubg=ubg,
        opts_ipopt=opts_ipopt,
        grad_f_func=grad_f_func,
        grad_g_func=grad_g_func,
        param_index=param_index,
        scale=scale,
        lower_band=lower_band,
        upper_band=upper_band,
        initial_guess=initial_guess,
        muscles=('ta', 'sol', 'gast'),
    )

    """
    grad_f_func = Function("grad_f", [w], [gradient(j, w)])
    grad_g_func = Function("grad_g", [w], [jacobian(g, w).T])

    online_callback = OnlineCallback(nx=w.shape[0],
                                     ng=g.shape[0],
                                     grad_f_func=grad_f_func,
                                     grad_g_func=grad_g_func)

    opts_ipopt = {
        "ipopt.max_iter": 2500,
        "ipopt.tol": 1e-4,
        "ipopt.print_info_string": "yes",
        "ipopt.linear_solver": "mumps", # ma57,
        "ipopt.hessian_approximation" : "exact",# ou "limited-memory"
        "iteration_callback": online_callback
    }



    opts_ipopt = {
        "ipopt.max_iter": 2500,
        "ipopt.tol": 1e-4,
        "ipopt.print_info_string": "yes",
        "ipopt.linear_solver": "mumps", # ma57,
        "ipopt.hessian_approximation" : "exact",# ou "limited-memory"
    }
    


    nlp = {'x': w, 'f': j, 'g': g}
    solver = nlpsol('solver', 'ipopt', nlp, opts_ipopt)

    sol = solver(x0=w0, lbx=lbw, ubx=ubw, lbg=lbg, ubg=ubg)
    """

    # ============ Extraction ============ #
    w_opt = sol['x'].full().flatten()
    cost = sol['f'].full().item()

    # [FIX 3] dé-scaler les n_opt premiers éléments
    param_opt = w_opt[:n_opt] * scale

    err_param = np.abs(muscle_tendon_parameters_num - param_opt)

    # ============ Reporting générique piloté par param_index ============ #
    _report_nlp_results(param_index, muscles=['ta', 'sol', 'gast'],
                        param_opt=param_opt,
                        ref=muscle_tendon_parameters_num,
                        err=err_param, cost=cost, n_trials=n_trials)

    return param_opt



def _report_nlp_results(param_index, muscles, param_opt, ref, err, cost,
                        n_trials):
    """Affichage des résultats, ordre et contenu dérivés de param_index."""
    # Unités lisibles selon le paramètre (affichage uniquement)
    UNIT = {
        'l0m':  ('m',   1.0),
        'lst':  ('m',   1.0),
        'phi0': ('deg', None),   # converti via rad2deg
        'f0m':  ('N',   1.0),
        'km':   ('-',   1.0),
        'kt':   ('-',   1.0),
    }

    def _fmt(p, idx):
        val_opt = param_opt[idx]
        val_ref = ref[idx]
        val_err = err[idx]
        unit, _ = UNIT.get(p, ('', 1.0))
        if p == 'phi0':
            return (f"{np.rad2deg(val_ref):8.2f} | {np.rad2deg(val_opt):8.2f} | "
                    f"{np.rad2deg(val_err):7.2f}  [deg]")
        return f"{val_ref:8.4f} | {val_opt:8.4f} | {val_err:7.4f}  [{unit}]"

    print("\n" + "=" * 70)
    print(f"Résultats NLP — n_trials = {n_trials} | coût = {cost:.4e}")
    print("=" * 70)
    print(f"{'Param':<10} {'Ref':>8} | {'Estimé':>8} | {'|Diff|':>7}")
    print("-" * 70)
    for p in param_index:                      # ordre = ordre d'optimisation
        for m in muscles:
            idx = param_index[p][m]
            print(f"{p+'_'+m:<10} {_fmt(p, idx)}")
    print("=" * 70)

"""
def optimization_nlp(data, initial_guess, lower_band, upper_band, skeleton_num,
                     muscle_tendon_parameters_num, unknown_parameters, casadi_function):


    Identify muscle-tendon parameters (ℓom, φo, Fom, ℓst) via an NLP.

    Args:
        data (np.ndarray): Shape (15, n_trials).
            Rows:
              [0]     : measured joint torque (N.m)
              [1:3]   : joint angles q (rad)
              [3:6]   : muscle activations [0, 1]
              [6:9]   : measured fiber lengths (m)
              [9:12]  : measured pennation angles (rad)
              [12:15] : measured tendon lengths (m)
        initial_guess, lower_band, upper_band : muscle-tendon parameters (size 12).
        skeleton_num : skeleton geometry (.osim).
        muscle_tendon_parameters_num : reference values used for comparison.
        unknown_parameters : SX vector of the 12 unknowns.
        casadi_function : dict of CasADi functions.

    Returns:
        xopt (np.ndarray): Shape (12,) — estimated (ℓom, φo, Fom, ℓst).


    # ============ Sanity checks ============ #
    assert data.shape[0] == 15, f"data must be (15, n_trials), got {data.shape}"

    initial_guess = np.asarray(initial_guess, dtype=float)
    lower_band = np.asarray(lower_band, dtype=float)
    upper_band = np.asarray(upper_band, dtype=float)


    assert np.all(lower_band <= upper_band), \
        f"lower_band > upper_band at: {np.where(lower_band > upper_band)[0]}"

    n_trials = data.shape[1]
    n_muscle = 3

    # ============ NLP set up ============ #
    # Decision variables, initial guess, and bounds
    w, w0, lbw, ubw = [], [], [], []
    j = 0
    # Constraints and their bounds
    g, lbg, ubg = [], [], []

    up = unknown_parameters

    # Per-trial residuals (kept for optional post-processing)
    e_torque, e_fiber, e_pennation = [], [], []

    # Unknown muscle-tendon parameters (12 values)
    w += [up]
    w0 += list(initial_guess)
    lbw += list(lower_band)
    ubw += list(upper_band)

    # Weights in the cost function
    w_torque = 1  # N.m
    w_length = 0.005  # mm
    w_angle = (1 / 180) * np.pi  # rad

    for trial in range(n_trials):
        # --- Extract data for this trial (one column = one trial) ---
        data_trial = data[:, trial]

        measured_torque = data_trial[0]  # N.m
        q_trial = [0, 0, 0, 0] + list(data_trial[1:3])  # rad
        a_trial = data_trial[3:6]  # [0, 1]
        measured_fiber_length = data_trial[6:9]  # m
        measured_pennation_angle = data_trial[9:12]  # rad
        measured_tendon_length = data_trial[12:15]  # m

        musculoskeletal_states_trial = q_trial + list(skeleton_num)
        neuromusculoskeletal_state_trial = np.concatenate(
            [a_trial, musculoskeletal_states_trial]
        )

        # MTU (muscle-tendon unit) length
        mtu_length = casadi_function['get_mtu_length'](musculoskeletal_states_trial)

        # --- Trial-specific decision variables ---
        tendon_length_k = SX.sym(f"Tendon_Length_{trial}", n_muscle)
        fiber_length_k = SX.sym(f"Fiber_Length_{trial}", n_muscle)
        pennation_angle_k = SX.sym(f"Pennation_Angle_{trial}", n_muscle)

        # --- Physiological bounds (robust to sign) ---
        fl_meas = np.abs(measured_fiber_length)
        tl_meas = np.abs(measured_tendon_length)
        pa_meas = np.abs(measured_pennation_angle)

        # Fiber length: positive, around the measured value
        lb_fl = fl_meas - fl_meas * 0.1
        ub_fl = fl_meas + fl_meas * 0.1

        # Pennation: physiological bounds, sign is free
        lb_pa = pa_meas - np.deg2rad(2)
        ub_pa = pa_meas + np.deg2rad(2)
        lb_pa = np.clip(lb_pa,0+epsilon,np.pi-epsilon)
        ub_pa = np.clip(ub_pa,0+epsilon,np.pi-epsilon)

        # Tendon length: positive, tight range
        lb_tl = tl_meas - tl_meas * 0.005
        ub_tl = tl_meas + tl_meas * 0.005

        lb_block = np.concatenate([lb_fl, lb_pa, lb_tl])
        ub_block = np.concatenate([ub_fl, ub_pa, ub_tl])

        # Initial guess, with lengths forced positive
        w0_k = np.concatenate([fl_meas, pa_meas, tl_meas])
        w_k = vertcat(fiber_length_k, pennation_angle_k, tendon_length_k)

        # Safety: check bounds consistency and that w0 lies within them
        assert np.all(lb_block <= ub_block), \
            f"Inverted bounds at trial {trial}"
        assert np.all((w0_k >= lb_block - 1e-9) & (w0_k <= ub_block + 1e-9)), \
            f"w0 out of bounds at trial {trial}: " \
            f"w0_k={w0_k}, lb={lb_block}, ub={ub_block}"

        w += [w_k]
        w0 += list(w0_k)
        lbw += list(lb_block)
        ubw += list(ub_block)

        # --- Equilibrium constraints (9 per trial) ---
        k = vertcat(a_trial, mtu_length, up)
        constraints = casadi_function['equilibrium_error_all_muscle'](w_k, k)
        g += [constraints]
        lbg += [0] * 9
        ubg += [0] * 9

        # --- Torque simulation ---
        all_states = vertcat(SX(neuromusculoskeletal_state_trial.tolist()), w_k)
        torque_simulated = casadi_function['get_joint_moment'](all_states, unknown_parameters)

        # --- Residuals ---
        e_torque_trials = measured_torque - torque_simulated
        e_fiber_trials = measured_fiber_length - fiber_length_k
        e_pennation_trials = measured_pennation_angle - pennation_angle_k

        j += w_torque * e_torque_trials ** 2
        j += sum1(w_length * e_fiber_trials ** 2)
        j += sum1(w_angle * e_pennation_trials ** 2)


        e_torque.append(e_torque_trials)
        e_fiber.append(e_fiber_trials)
        e_pennation.append(e_pennation_trials)

    # ============ Assembly ============ #
    w = vertcat(*w)
    g = vertcat(*g)
    print("J shape:", j.shape)
    print(f"n_trials : {n_trials} | size(w) : {w.shape} | size(g) : {g.shape}")

    w0 = np.array(w0, dtype=float)
    lbw = np.array(lbw, dtype=float)
    ubw = np.array(ubw, dtype=float)
    lbg = np.array(lbg, dtype=float)
    ubg = np.array(ubg, dtype=float)

    # Global checks
    if np.any(np.isnan(w0)):
        print('NaNs in w0 at:', np.where(np.isnan(w0)))
        return None
    if np.any(np.isinf(w0)):
        print('Infs in w0 at:', np.where(np.isinf(w0)))
        return None
    bad = np.where(lbw > ubw)[0]
    if len(bad):
        print(f'lbw > ubw at indices: {bad[:20]}...')
        return None

    print('w0 is valid')

    # ============ NLP solver ============ #

    opts_ipopt = {
        "ipopt.max_iter": 2500
    }

    nlp = {'x': w, 'f': j, 'g': g}
    solver = nlpsol('solver', 'ipopt', nlp,opts_ipopt)

    print(solver)

    sol = solver(x0=w0, lbx=lbw, ubx=ubw, lbg=lbg, ubg=ubg)

    # ============ Extraction ============ #
    w_opt = sol['x'].full().flatten()
    cost = sol['f'].full().item()

    param_opt = w_opt[:12]
    err_param = np.abs(muscle_tendon_parameters_num - param_opt)

    print("\n\nDifference between input and estimated muscle-tendon parameters:")
    print(f"Number of trials : {n_trials}")
    print(f"diff ℓom : {err_param[0:3]}")
    print(f"diff φo  : {err_param[3:6]}")
    print(f"diff Fom : {err_param[6:9]}")
    print(f"diff ℓst : {err_param[9:12]}")

    print("input muscle-tendon parameters:")
    print(f"ℓom : {muscle_tendon_parameters_num[0:3]}")
    print(f"φo  : {muscle_tendon_parameters_num[3:6]}")
    print(f"Fom : {muscle_tendon_parameters_num[6:9]}")
    print(f"ℓst : {muscle_tendon_parameters_num[9:12]}")

    print("Estimated muscle-tendon parameters:")
    print(f"Cost : {cost}")
    print(f"ℓom : {param_opt[0:3]}")
    print(f"φo  : {param_opt[3:6]}")
    print(f"Fom : {param_opt[6:9]}")
    print(f"ℓst : {param_opt[9:12]}")


    return param_opt


"""


def generate_estimated_data(data, skeleton_num, muscle_tendon_parameters_num,
                            casadi_function,
                            output_dir='data_estime', filename='data_estime.csv',
                            save_npy=True, save_csv=True, verbose=True):
    """
    Génère un dataset "estimé" par le modèle direct à partir d'un dataset d'entrée.

    Pour chaque trial de `data`, on réutilise les angles articulaires et les
    activations mesurées, puis on résout l'équilibre muscle-tendon pour
    recalculer :
        - les longueurs de fibres
        - les angles de pennation
        - les longueurs de tendon
        - le torque articulaire
    avec les paramètres muscle-tendon fournis.

    Utile pour :
        - valider la calibration (mesure vs modèle sur mêmes q, a)
        - générer un jeu de référence pour évaluer l'optim

    Args:
        data (np.ndarray): Shape (15, n_trials). Seules les lignes [1:3] (q)
            et [3:6] (a) sont utilisées ; les autres sont recalculées.
            Unités : q en rad, a ∈ [0, 1].
        skeleton_num : géométrie squelette (.osim scalé).
        muscle_tendon_parameters_num (array-like, 12): (ℓom, φo, Fom, ℓst) × 3 muscles.
            Ordre : [TA, SOL, GAST] pour chaque paramètre.
        casadi_function (dict): fonctions CasADi.
        output_dir (str): dossier de sortie. Créé si inexistant.
        filename (str): nom de base du fichier.
        save_npy, save_csv (bool): formats de sauvegarde.
        verbose (bool): logs de progression.

    Returns:
        header (list[str]): noms des lignes.
        data_est (np.ndarray): Shape (15, n_success_trials), mêmes unités que `data`
            (angles en rad, longueurs en m, torque en N.m).
        path_csv (str | None): chemin du CSV écrit.
    """

    header = [
        'ankle_torque', 'q_knee', 'q_ankle',
        'a_tibialis', 'a_soleus', 'a_gastrocnemius',
        'fiber_length_tibialis', 'fiber_length_soleus', 'fiber_length_gastrocnemius',
        'pennation_angle_tibialis', 'pennation_angle_soleus', 'pennation_angle_gastrocnemius',
        'tendon_length_tibialis', 'tendon_length_soleus', 'tendon_length_gastrocnemius'
    ]

    # --- Validation input ---
    assert data.shape[0] == 15, f"data doit être (15, n_trials), reçu {data.shape}"
    n_trials = data.shape[1]

    # Paramètres muscle par muscle
    mtp = np.array(muscle_tendon_parameters_num)
    mtp_ta = mtp[[0, 3, 6, 9]]
    mtp_sol = mtp[[1, 4, 7, 10]]
    mtp_gast = mtp[[2, 5, 8, 11]]

    data_est = np.zeros((15, n_trials))
    n_fail, n_ok = 0, 0

    if verbose:
        print(f'\n================== Generating estimated data ==================')
        print(f"  n_trials input : {n_trials}")

    for trial in range(n_trials):
        if verbose and (trial + 1) % max(1, n_trials // 20) == 0:
            print(f"  progress : {100 * (trial + 1) / n_trials:.0f} %")

        col = data[:, trial]

        # --- Récupération des inputs (q, a) ---
        q_trial = col[1:3]  # [q_knee, q_ankle] en rad
        a_trial = col[3:6]  # [a_TA, a_SOL, a_GAST]
        a_ta, a_sol, a_gast = a_trial

        # État musculosquelettique complet
        q_num = [0, 0, 0, 0] + list(q_trial)
        musc_states = q_num + list(skeleton_num)
        neuromusc_state = np.concatenate([a_trial, musc_states])

        # Longueur MTU
        mtu_length = casadi_function['get_mtu_length'](musc_states)
        l_ta = float(mtu_length[0])
        l_sol = float(mtu_length[1])
        l_gast = float(mtu_length[2])

        # --- Équilibre muscle-tendon par muscle ---
        x_ta, s_ta, cost = root_muscle_dynamics(a_ta, l_ta, mtp_ta, 'tibialis', casadi_function)
        x_sol, s_sol, cost = root_muscle_dynamics(a_sol, l_sol, mtp_sol, 'soleus', casadi_function)
        x_gast, s_gast, cost = root_muscle_dynamics(a_gast, l_gast, mtp_gast, 'gastrocnemius', casadi_function)

        """
        if 'fail' in (s_ta, s_sol, s_gast):
            n_fail += 1
            continue
            
        """

        # Variables rootées
        fl = np.array([x_ta[0, 0], x_sol[0, 0], x_gast[0, 0]]).flatten()
        penn = np.array([x_ta[1, 0], x_sol[1, 0], x_gast[1, 0]]).flatten()
        tl = np.array([x_ta[2, 0], x_sol[2, 0], x_gast[2, 0]]).flatten()

        rooted = np.concatenate([fl, penn, tl])
        all_state = np.concatenate([neuromusc_state, rooted])

        # Torque estimé
        ankle_torque = float(casadi_function['get_joint_moment'](all_state, mtp))

        # --- Remplissage (mêmes unités que `data`) ---
        data_est[:, n_ok] = [
            ankle_torque,
            q_trial[0], q_trial[1],  # rad
            a_ta, a_sol, a_gast,
            fl[0], fl[1], fl[2],  # m
            penn[0], penn[1], penn[2],  # rad
            tl[0], tl[1], tl[2],  # m
        ]
        n_ok += 1

    data_est = data_est[:, :n_ok]

    if verbose:
        print(f"  → {n_ok} / {n_trials} trials OK ({100 * n_fail / n_trials:.1f} % fail)")

    # --- Sauvegarde ---
    os.makedirs(output_dir, exist_ok=True)
    base = os.path.splitext(filename)[0]
    path_csv = None

    if save_csv:
        # Pour le CSV : transposer en (n_trials, 15) + angles en degrés pour lisibilité
        path_csv = os.path.join(output_dir, base + '.csv')
        df = pd.DataFrame(data_est.T, columns=header)
        df['q_knee'] = np.rad2deg(df['q_knee'])
        df['q_ankle'] = np.rad2deg(df['q_ankle'])
        for m in ['tibialis', 'soleus', 'gastrocnemius']:
            df[f'pennation_angle_{m}'] = np.rad2deg(df[f'pennation_angle_{m}'])
        df.to_csv(path_csv, index=False)
        if verbose:
            print(f"  → CSV (deg, n_trials × 15) : {path_csv}")

    if save_npy:
        # NPY : shape (15, n_trials) en unités SI (rad, m, N.m) — directement utilisable
        path_npy = os.path.join(output_dir, base + '.npy')
        np.save(path_npy, data_est)
        if verbose:
            print(f"  → NPY (rad, 15 × n_trials) : {path_npy}")

    return header, data_est, path_csv

def save_data_to_xlsx(data, header_data, folder, name,
                      use_symbolic=True, convert_units=True):
    """
    Save data matrix to Excel with custom column order and symbolic headers.

    Args:
        data : np.ndarray, shape (n_features, n_trials)
        header_data : list of str — names of the rows of data.
        folder, name : output directory and filename (no extension).
        use_symbolic : if True, Excel headers use symbols (τ_ankle, ...).
        convert_units : if True, converts m -> cm and rad -> deg.

    Returns:
        str : full path to the saved Excel file.
    """

    # Mapping : explicit name -> symbol (for Excel headers)
    SYMBOL_MAP = {
        'ankle_torque': 'τ_ankle',
        'q_knee': 'q_1',
        'q_ankle': 'q_2',
        'a_tibialis': 'a_Ta',
        'a_soleus': 'a_sol',
        'a_gastrocnemius': 'a_Gast',
        'fiber_length_tibialis': 'lf_Ta^',
        'fiber_length_soleus': 'lf_sol^',
        'fiber_length_gastrocnemius': 'lf_Gast^',
        'pennation_angle_tibialis': 'φ_Ta',
        'pennation_angle_soleus': 'φ_sol',
        'pennation_angle_gastrocnemius': 'φ_Gast',
        'tendon_length_tibialis': 'lt_Ta',
        'tendon_length_soleus': 'lt_sol',
        'tendon_length_gastrocnemius': 'lt_Gast',
    }

    # Desired column order in Excel (grouped per muscle)
    output_order = [
        'ankle_torque',
        'q_knee', 'q_ankle',
        'a_tibialis',      'fiber_length_tibialis',       'pennation_angle_tibialis',
        'a_soleus',        'fiber_length_soleus',         'pennation_angle_soleus',
        'a_gastrocnemius', 'fiber_length_gastrocnemius',  'pennation_angle_gastrocnemius',
    ]

    ANGLE_ROWS_RAD_TO_DEG = [
        'q_knee', 'q_ankle',
        'pennation_angle_tibialis',
        'pennation_angle_soleus',
        'pennation_angle_gastrocnemius',
    ]
    LENGTH_ROWS_M_TO_CM = [
        'fiber_length_tibialis',
        'fiber_length_soleus',
        'fiber_length_gastrocnemius',
        'tendon_length_tibialis',
        'tendon_length_soleus',
        'tendon_length_gastrocnemius',
    ]

    # --- Validations --- #
    if len(header_data) != data.shape[0]:
        raise ValueError(
            f"Length mismatch: header_data has {len(header_data)} entries "
            f"but data has {data.shape[0]} rows."
        )

    missing = [n for n in output_order if n not in header_data]
    if missing:
        raise ValueError(
            f"output_order contains names not in header_data: {missing}"
        )

    # --- Réordonnancement des lignes selon output_order --- #
    row_indices = [header_data.index(n) for n in output_order]
    data_reordered = data[row_indices, :].astype(float).copy()

    # --- Conversion d'unités --- #
    if convert_units:
        for i, varname in enumerate(output_order):
            if varname in ANGLE_ROWS_RAD_TO_DEG:
                data_reordered[i, :] = np.rad2deg(data_reordered[i, :])
            elif varname in LENGTH_ROWS_M_TO_CM:
                data_reordered[i, :] = data_reordered[i, :] * 100.0

    # --- Labels Excel (une seule fois) --- #
    columns = [SYMBOL_MAP[n] for n in output_order] if use_symbolic else list(output_order)

    # --- Sauvegarde --- #
    df = pd.DataFrame(data_reordered.T, columns=columns)
    os.makedirs(folder, exist_ok=True)
    excel_path = os.path.join(folder, f"{name}.xlsx")
    df.to_excel(excel_path, index=False)

    print(f'data saved in : {excel_path}')
    return excel_path

def load_data_from_xlsx(folder, name, header_data,
                        use_symbolic=True, convert_units=True):
    """
    Load data matrix from Excel, reversing the transformations applied
    by save_data_to_xlsx.

    Variables listed in `header_data` but absent from the Excel file are
    filled with NaN values, preserving the matrix shape (n_features, n_trials).

    Args:
        folder : str — directory containing the file.
        name : str — base filename (no extension).
        header_data : list of str — explicit names of the rows requested.
        use_symbolic : bool, default True
            If True, expects Excel columns to use symbolic names.
        convert_units : bool, default True
            If True, reverses cm -> m and deg -> rad conversion.

    Returns:
        data : np.ndarray, shape (len(header_data), n_trials)
            Matrix with rows in the order of header_data. Missing variables
            are filled with NaN.
    """

    SYMBOL_MAP = {
        'ankle_torque': 'τ_ankle',
        'q_knee': 'q_1',
        'q_ankle': 'q_2',
        'a_tibialis': 'a_Ta',
        'a_soleus': 'a_sol',
        'a_gastrocnemius': 'a_Gast',
        'fiber_length_tibialis': 'lf_Ta^',
        'fiber_length_soleus': 'lf_sol^',
        'fiber_length_gastrocnemius': 'lf_Gast^',
        'pennation_angle_tibialis': 'φ_Ta',
        'pennation_angle_soleus': 'φ_sol',
        'pennation_angle_gastrocnemius': 'φ_Gast',
        'tendon_length_tibialis': 'lt_Ta',
        'tendon_length_soleus': 'lt_sol',
        'tendon_length_gastrocnemius': 'lt_Gast',
    }

    ANGLE_ROWS_RAD_TO_DEG = [
        'q_knee', 'q_ankle',
        'pennation_angle_tibialis',
        'pennation_angle_soleus',
        'pennation_angle_gastrocnemius',
    ]
    LENGTH_ROWS_M_TO_CM = [
        'fiber_length_tibialis',
        'fiber_length_soleus',
        'fiber_length_gastrocnemius',
        'tendon_length_tibialis',
        'tendon_length_soleus',
        'tendon_length_gastrocnemius',
    ]

    # --- Lecture du fichier Excel --- #
    if name.endswith('.xlsx'):
        excel_path = os.path.join(folder, name)
    else:
        excel_path = os.path.join(folder, f"{name}.xlsx")

    if not os.path.exists(excel_path):
        raise FileNotFoundError(f"Excel file not found: {excel_path}")
    df = pd.read_excel(excel_path)

    # --- Reverse mapping symbole -> nom explicite --- #
    if use_symbolic:
        reverse_symbol_map = {sym: name for name, sym in SYMBOL_MAP.items()}
        excel_to_explicit = {col: reverse_symbol_map.get(col, col) for col in df.columns}
    else:
        excel_to_explicit = {col: col for col in df.columns}

    explicit_to_excel = {v: k for k, v in excel_to_explicit.items()}
    available_explicit = set(explicit_to_excel.keys())

    # --- Construction de la matrice avec NaN pour les variables manquantes --- #
    n_trials = len(df)
    n_features = len(header_data)
    data = np.full((n_features, n_trials), np.nan)

    missing = []
    for i, varname in enumerate(header_data):
        if varname in available_explicit:
            excel_col = explicit_to_excel[varname]
            data[i, :] = df[excel_col].to_numpy(dtype=float)
        else:
            missing.append(varname)
        # sinon : ligne reste à NaN

    # --- Conversion inverse d'unités (uniquement sur les lignes non-NaN) --- #
    if convert_units:
        for i, varname in enumerate(header_data):
            if varname in missing:
                continue  # ligne NaN, pas de conversion
            if varname in ANGLE_ROWS_RAD_TO_DEG:
                data[i, :] = np.deg2rad(data[i, :])
            elif varname in LENGTH_ROWS_M_TO_CM:
                data[i, :] = data[i, :] / 100.0

    if missing:
        print(f"[load_data_from_xlsx] Filled with NaN (not in file): {missing}")
    print(f'data loaded from : {excel_path} '
          f'({n_features - len(missing)}/{n_features} variables)')

    return data

def compare_datasets(data_a, data_b, header_data, label_a='A', label_b='B'):
    """
    Compare two data matrices row-by-row and report statistics on the
    differences (per variable).

    Parameters
    ----------
    data_a, data_b : np.ndarray, shape (n_features, n_trials)
        The two datasets to compare. Must have identical shape.
    header_data : list of str
        Variable names corresponding to the rows of both matrices.
    label_a, label_b : str, optional
        Names used in the printed report (e.g. 'original', 'reloaded').

    Returns
    -------
    pd.DataFrame
        Per-variable statistics with columns :
            mean_diff, std_diff, min_diff, max_diff,
            mean_abs_diff, max_abs_diff, rmse, n_trials.
    """
    # --- Validations --- #
    if data_a.shape != data_b.shape:
        raise ValueError(
            f"Shape mismatch: {label_a} has {data_a.shape}, "
            f"{label_b} has {data_b.shape}."
        )
    if len(header_data) != data_a.shape[0]:
        raise ValueError(
            f"header_data has {len(header_data)} entries but matrices have "
            f"{data_a.shape[0]} rows."
        )

    # --- Calcul des différences --- #
    diff = data_a.astype(float) - data_b.astype(float)
    abs_diff = np.abs(diff)

    plt.imshow(diff, aspect='auto', cmap='hot')
    plt.colorbar(label='|diff|')
    plt.yticks(range(len(header_data)), header_data)
    plt.xlabel('Trial')
    plt.title('Round-trip absolute difference')
    plt.tight_layout()
    plt.show()

    stats = pd.DataFrame({
        'mean_diff':     diff.mean(axis=1),
        'std_diff':      diff.std(axis=1, ddof=1) if diff.shape[1] > 1 else np.zeros(diff.shape[0]),
        'min_diff':      diff.min(axis=1),
        'max_diff':      diff.max(axis=1),
        'mean_abs_diff': abs_diff.mean(axis=1),
        'max_abs_diff':  abs_diff.max(axis=1),
        'rmse':          np.sqrt((diff ** 2).mean(axis=1)),
        'n_trials':      diff.shape[1],
    }, index=header_data)

    # --- Rapport console --- #
    print(f"\n{'='*70}")
    print(f"  Comparison : {label_a} vs {label_b}")
    print(f"  Shape : {data_a.shape}  ({data_a.shape[0]} variables × {data_a.shape[1]} trials)")
    print(f"{'='*70}")

    # Format colonné lisible
    with pd.option_context('display.float_format', '{:.3e}'.format,
                           'display.max_rows', None,
                           'display.width', None):

        print(stats)

    # Résumé global
    overall_max = abs_diff.max()
    overall_mean = abs_diff.mean()
    print(f"\nOverall max |diff| : {overall_max:.3e}")
    print(f"Overall mean |diff|: {overall_mean:.3e}")

    # Détection de divergence : variable avec la plus grande erreur relative
    nonzero_mask = np.abs(data_a).max(axis=1) > 1e-12
    if nonzero_mask.any():
        rel_err = np.where(
            nonzero_mask,
            abs_diff.max(axis=1) / np.maximum(np.abs(data_a).max(axis=1), 1e-12),
            0.0,
        )
        worst_idx = int(np.argmax(rel_err))
        print(f"\nWorst relative error : '{header_data[worst_idx]}' "
              f"({rel_err[worst_idx]:.3e})")

    print(f"{'='*70}\n")
    return stats

def get_initial_guess(muscle_tendon_parameters_num, data, param_config,
                      strategy='measured', verbose=True):
    """
    Génère initial_guess, lower_band et upper_band pour la calibration
    muscle-tendon, uniquement pour les paramètres déclarés 'sym' dans
    param_config.

    Conventions
    -----------
    param_config : dict ordonné {nom_param: 'sym' | 'fixed'} parmi
        'l0m', 'phi0', 'f0m', 'lst', 'km', 'kt'.
        - 'sym'   : paramètre optimisé → inclus dans le vecteur de sortie.
        - 'fixed' : paramètre figé → exclu du vecteur de sortie, mais sa
                    valeur reste présente dans muscle_tendon_parameters_num.
        L'ordre des clés détermine à la fois :
          (1) le layout de muscle_tendon_parameters_num (array de référence) :
              un bloc [ta, sol, gast] par clé, dans l'ordre des clés ;
          (2) l'ordre des paramètres 'sym' dans le vecteur optimisé.

    muscle_tendon_parameters_num : np.ndarray plat, longueur 3*len(param_config).
        Valeurs de référence scalées (modèle OpenSim), bloc [ta, sol, gast]
        par paramètre, dans l'ordre des clés de param_config.

    data : np.ndarray, shape (15, ntrials). Mesures écho :
        indices 6-8  = fascicle length (ta, sol, gast)
        indices 9-11 = pennation angle (ta, sol, gast)
        indices 12-14= tendon length   (ta, sol, gast)

    Stratégies
    ----------
    - 'literature' : bornes centrées sur les valeurs de littérature (Rajagopal 2015).
    - 'scaled'     : bornes centrées sur muscle_tendon_parameters_num.
    - 'hybrid'     : init = scaled, bornes = scaled avec garde-fous littérature.
    - 'measured'   : règles data-driven explicites :
                       l0m  : init = mean(FL),  bornes = [min(FL), max(FL)]
                       phi0 : init = mean(PA),  bornes = [min(PA), max(PA)]
                       f0m  : init = num,        bornes = [0.5*num, 2.0*num]
                       lst  : init = min(TL),    bornes = [min(TL)-0.1, min(TL)+0.02]
                       km   : init = num,        bornes = [0.5*num, 2.0*num]
                       kt   : init = num,        bornes = [0.5*num, 2.0*num]

    Returns
    -------
    initial_guess, upper_band, lower_band : np.ndarray, shape (3*n_sym,)
    param_index : dict {nom_param: {muscle: idx}} pour relire le vecteur optimisé.
    """

    muscles = ['ta', 'sol', 'gast']
    ALL_PARAMS = ['l0m', 'phi0', 'f0m', 'lst', 'km', 'kt']

    # === Indices des MESURES (data) ===
    IDX_FL = {'ta': 6,  'sol': 7,  'gast': 8}   # fascicle length
    IDX_PA = {'ta': 9,  'sol': 10, 'gast': 11}  # pennation angle
    IDX_TL = {'ta': 12, 'sol': 13, 'gast': 14}  # tendon length

    # === Valeurs de référence (Rajagopal 2015) ===
    LITERATURE = {
        'l0m':  {'ta': 0.0683, 'sol': 0.0440, 'gast': 0.0600},
        'phi0': {'ta': np.deg2rad(9.6),
                 'sol': np.deg2rad(28.3),
                 'gast': np.deg2rad(9.9)},
        'lst':  {'ta': 0.2230, 'sol': 0.2450, 'gast': 0.3800},
    }

    # === Marges (literature / scaled / hybrid) ===
    MARGINS = {
        'l0m':  {'low': 0.70, 'high': 1.30},
        'phi0': {'low_rad': np.deg2rad(-10), 'high_rad': np.deg2rad(10)},
        'f0m':  {'low': 0.10, 'high': 2.5},
        'lst':  {'low': 0.80, 'high': 1.25},
    }

    # === Facteurs scalés f0m / km / kt ===
    SCALED_LOW  = 0.5
    SCALED_HIGH = 3.0

    # === Marges lst en mètres (measured) ===
    LST_MARGIN_LOW  = 0.20
    LST_MARGIN_HIGH = 0.03

    # ------------------------------------------------------------------
    # Validation de param_config et de la stratégie
    # ------------------------------------------------------------------
    if not isinstance(param_config, dict):
        raise TypeError("param_config doit être un dict {nom: 'sym'|'fixed'}.")

    if strategy not in ('literature', 'scaled', 'hybrid', 'measured'):
        raise ValueError(f"Strategy inconnue : {strategy}")

    ref_params = list(param_config.keys())          # ordre du modèle / array
    unknown = [p for p in ref_params if p not in ALL_PARAMS]
    if unknown:
        raise ValueError(f"Paramètres inconnus dans param_config : {unknown}")

    active = [p for p in ref_params if param_config[p] == 'sym']  # à optimiser

    if not active:
        raise ValueError("Aucun paramètre 'sym' dans param_config.")

    # ------------------------------------------------------------------
    # Indices de l'array de référence, dérivés des paramètres 'sym'
    # (l'array ne contient QUE les 'sym', un bloc [ta, sol, gast] chacun,
    #  dans l'ordre d'apparition dans param_config)
    # ------------------------------------------------------------------
    expected_len = 3 * len(active)
    if len(muscle_tendon_parameters_num) != expected_len:
        raise ValueError(
            f"muscle_tendon_parameters_num a {len(muscle_tendon_parameters_num)} "
            f"éléments mais {len(active)} paramètres 'sym' ({active}) en attendent "
            f"{expected_len} (un bloc [ta,sol,gast] par paramètre 'sym').")

    REF_IDX = {}
    j = 0
    for p in active:
        REF_IDX[p] = {'ta': j, 'sol': j + 1, 'gast': j + 2}
        j += 3

    # ------------------------------------------------------------------
    # Indices du vecteur optimisé (paramètres 'sym' uniquement, dans l'ordre)
    # ------------------------------------------------------------------
    param_index = {}
    idx = 0
    for p in active:
        param_index[p] = {}
        for m in muscles:
            param_index[p][m] = idx
            idx += 1
    n_params = idx

    initial_guess = np.zeros(n_params)
    lower_band    = np.zeros(n_params)
    upper_band    = np.zeros(n_params)

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------
    def _clean(arr):
        arr = arr[~np.isnan(arr)]
        if arr.size == 0:
            raise ValueError("Aucune mesure valide (que des NaN).")
        return arr

    def _num(p, m):
        """Valeur de référence scalée pour le paramètre p, muscle m."""
        if p not in REF_IDX:
            raise KeyError(
                f"'{p}' est 'sym' mais absent de muscle_tendon_parameters_num.")
        return muscle_tendon_parameters_num[REF_IDX[p][m]]

    def _compute_measured(p, m):
        """Règles data-driven explicites de la stratégie 'measured'."""
        if p == 'l0m':
            fl = _clean(data[IDX_FL[m], :])
            return np.mean(fl), np.min(fl), np.max(fl)
        if p == 'phi0':
            pa = _clean(data[IDX_PA[m], :])
            return np.mean(pa), np.min(pa), np.max(pa)
        if p == 'lst':
            tl = _clean(data[IDX_TL[m], :])
            tl_min = np.min(tl)
            return tl_min, tl_min - LST_MARGIN_LOW, tl_min + LST_MARGIN_HIGH
        # f0m, km, kt : scalé ±facteur
        v = _num(p, m)
        return v, SCALED_LOW * v, SCALED_HIGH * v

    def _compute_centered(p, m):
        """Stratégies literature / scaled / hybrid."""
        if strategy == 'literature' and p in LITERATURE:
            center = LITERATURE[p][m]
        else:  # scaled, hybrid, ou paramètre sans référence littérature
            center = _num(p, m)

        if p == 'l0m':
            lo = center * MARGINS['l0m']['low']
            hi = center * MARGINS['l0m']['high']
            if strategy == 'hybrid':
                lo = max(lo, LITERATURE['l0m'][m] * 0.60)
                hi = min(hi, LITERATURE['l0m'][m] * 1.40)
            return center, lo, hi
        if p == 'phi0':
            lo = max(center + MARGINS['phi0']['low_rad'], np.deg2rad(2))
            hi = min(center + MARGINS['phi0']['high_rad'], np.deg2rad(40))
            return center, lo, hi
        if p == 'f0m':
            return (center,
                    center * MARGINS['f0m']['low'],
                    center * MARGINS['f0m']['high'])
        if p == 'lst':
            return (center,
                    center * MARGINS['lst']['low'],
                    center * MARGINS['lst']['high'])
        # km, kt : scalé ±facteur (pas de logique littérature)
        v = _num(p, m)
        return v, SCALED_LOW * v, SCALED_HIGH * v

    # ------------------------------------------------------------------
    # Remplissage du vecteur optimisé
    # ------------------------------------------------------------------
    for p in active:
        for m in muscles:
            if strategy == 'measured':
                init, lo, hi = _compute_measured(p, m)
            else:
                init, lo, hi = _compute_centered(p, m)
            i = param_index[p][m]
            initial_guess[i] = init
            lower_band[i]    = lo
            upper_band[i]    = hi

    # ------------------------------------------------------------------
    # Diagnostic
    # ------------------------------------------------------------------
    if verbose:
        _print_initial_guess(initial_guess, lower_band, upper_band,
                             muscles, param_index, active, strategy)
        # cohérence géométrique seulement si l0m, phi0 et lst sont tous 'sym'
        if all(p in param_index for p in ('l0m', 'phi0', 'lst')):
            _check_geometric_consistency(
                data, IDX_FL, IDX_PA, IDX_TL, initial_guess,
                param_index['l0m'], param_index['phi0'], param_index['lst'],
                muscles)

    return initial_guess, upper_band, lower_band, param_index


def _print_initial_guess(init, lo, hi, muscles, param_index, active, strategy):
    """Affichage formaté de l'initial guess et des bornes."""
    print("=" * 78)
    print(f"Initial guess — stratégie '{strategy}' — params : {', '.join(active)}")
    print("=" * 78)
    print(f"{'Param':<14} {'Init':>10} {'Lower':>10} {'Upper':>10} "
          f"{'Marge L':>9} {'Marge U':>9}")
    print("-" * 78)

    for p in active:
        for m in muscles:
            i = param_index[p][m]
            label = f"{p}_{m}"
            margin_lo = (init[i] - lo[i]) / init[i] * 100 if init[i] != 0 else 0
            margin_hi = (hi[i] - init[i]) / init[i] * 100 if init[i] != 0 else 0
            print(f"{label:<14} {init[i]:>10.4f} {lo[i]:>10.4f} {hi[i]:>10.4f} "
                  f"{margin_lo:>8.1f}% {margin_hi:>8.1f}%")
    print("=" * 78)


def _check_geometric_consistency(data, IDX_FL, IDX_PA, IDX_TL,
                                 init, IDX_LOM, IDX_PHI, IDX_LST, muscles):
    """
    Vérifie que les mesures (ℓ^M, α, ℓ^T) sont cohérentes avec
    les valeurs initiales des paramètres et avec une longueur muscle-tendon
    physiologiquement plausible.
    """
    print("\nDiagnostic de cohérence des mesures :")
    print("-" * 78)
    warnings = []

    for m in muscles:
        fl = data[IDX_FL[m], :]
        pa = data[IDX_PA[m], :]
        tl = data[IDX_TL[m], :]
        lmt_reconstructed = fl * np.cos(pa) + tl

        fl_mean, fl_std = np.nanmean(fl), np.nanstd(fl)
        pa_mean, pa_std = np.nanmean(pa), np.nanstd(pa)
        tl_mean, tl_std = np.nanmean(tl), np.nanstd(tl)
        lmt_mean = np.nanmean(lmt_reconstructed)
        lmt_range = np.nanmax(lmt_reconstructed) - np.nanmin(lmt_reconstructed)

        lom_init = init[IDX_LOM[m]]
        lst_init = init[IDX_LST[m]]
        lmt_expected = lom_init * np.cos(init[IDX_PHI[m]]) + lst_init

        print(f"\n  [{m.upper()}]")
        print(f"    Mesures  : ℓ^M = {fl_mean*100:.2f}±{fl_std*100:.2f} cm   "
              f"α = {np.rad2deg(pa_mean):.1f}±{np.rad2deg(pa_std):.1f}°   "
              f"ℓ^T = {tl_mean*100:.2f}±{tl_std*100:.2f} cm")
        print(f"    ℓ^MT reconstruit (mesures)      : {lmt_mean*100:.2f} cm "
              f"(amplitude {lmt_range*100:.2f} cm)")
        print(f"    ℓ^MT attendu (params init)      : {lmt_expected*100:.2f} cm")
        print(f"    Δ (mesures - params)            : "
              f"{(lmt_mean - lmt_expected)*100:+.2f} cm")

        # Alertes
        if abs(lmt_mean - lmt_expected) > 0.02:  # >2 cm de désaccord
            warnings.append(f"  ⚠ {m.upper()} : ℓ^MT mesuré et ℓ^MT issu des params "
                          f"diffèrent de {(lmt_mean - lmt_expected)*100:+.1f} cm. "
                          f"Vérifier le modèle géométrique de longueur muscle-tendon.")

        # Le tendon ne devrait jamais être plus court que lst
        if np.nanmin(tl) < lst_init * 0.95:
            warnings.append(f"  ⚠ {m.upper()} : ℓ^T mesuré "
                          f"({np.nanmin(tl)*100:.1f} cm) < lst_init "
                          f"({lst_init*100:.1f} cm). Tendon mesuré trop court.")

        # La fiber length devrait osciller autour de lom (plage typique 0.5-1.5 ℓ_0^M)
        if fl_mean < lom_init * 0.5 or fl_mean > lom_init * 1.5:
            warnings.append(f"  ⚠ {m.upper()} : ℓ^M moyen "
                          f"({fl_mean*100:.1f} cm) très éloigné de lom_init "
                          f"({lom_init*100:.1f} cm).")

    if warnings:
        print("\nAlertes :")
        for w in warnings:
            print(w)
        print("\n→ Si Δ > 2 cm, l'optimisation NLP risque l'infaisabilité.")
        print("  Vérifier en priorité le calcul de ℓ^MT(q_knee, q_ankle).")
    else:
        print("\n  ✓ Pas d'incohérence majeure détectée.")
    print("=" * 78)
