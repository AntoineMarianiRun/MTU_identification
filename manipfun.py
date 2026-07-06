# function to work on the data from neuro-musculo-skeletal (emg and echo driven)
#################################################################################
# lib
import pandas as pd
import numpy as np
import tkinter as tk
from tkinter import filedialog, messagebox
from pathlib import Path
import matplotlib.pyplot as plt
import os
from datetime import datetime

# function
def select_folder_and_get_files():
    root = tk.Tk()
    root.withdraw()

    folder = filedialog.askdirectory(
        title="select the interest folder"
    )

    if not folder:
        return None, None, None

    folder_path = Path(folder)
    results = {
        'osim': None,
        'train': None,
        'test': None
    }

    # find .osim
    osim_files = list(folder_path.glob('*.osim'))
    if osim_files:
        results['osim'] = osim_files[0]

    # Chercher Excel (train et test)
    excel_files = list(folder_path.glob('*.xlsx')) + \
                  list(folder_path.glob('*.xls')) + \
                  list(folder_path.glob('*.csv'))

    for excel_file in excel_files:
        filename_lower = excel_file.name.lower()
        if 'train' in filename_lower and not results['train']:
            results['train'] = excel_file
        elif 'test' in filename_lower and not results['test']:
            results['test'] = excel_file

    # return paht as string
    return (str(results['osim']) if results['osim'] else None,
            str(results['train']) if results['train'] else None,
            str(results['test']) if results['test'] else None)

def split_path_name(full_path):
    """
    Avoir le nom et le path
        Args:
        full_path (str):
            [full_path]
    """
    full_path = Path(full_path)
    folder_name = str(full_path.parent)
    file_name = full_path.name
    return folder_name, file_name

def import_data_from_excel(full_path):
    """
    Importe les données depuis un fichier Excel

    Args:
        full_path (str): Chemin complet du fichier Excel

    Returns:
        data (np.ndarray): Shape (15, ntrials) contenant :
                          [ankle_torque, q_knee, q_ankle,
                           a_tibialis, a_soleus, a_gastrocnemius,
                           fiber_length_tibialis, fiber_length_soleus, fiber_length_gastrocnemius,
                           pennation_angle_tibialis, pennation_angle_soleus, pennation_angle_gastrocnemius,
                           tendon_length_tibialis, tendon_length_soleus, tendon_length_gastrocnemius]
    """

    # Vérifier que le fichier existe
    if not Path(full_path).exists():
        raise FileNotFoundError(f"Le fichier {full_path} n'existe pas")

    # Charger le fichier
    print(f"Chargement : {full_path}")
    try:
        if full_path.endswith('.csv'):
            df = pd.read_csv(full_path)
        else:
            df = pd.read_excel(full_path)
    except Exception as e:
        raise ValueError(f"Erreur lecture fichier : {e}")

    # Mapping des colonnes Excel vers les noms attendus
    column_mapping = {
        'τ_ankle': 'ankle_torque',
        'q_1': 'q_knee',
        'q_2': 'q_ankle',
        'a_Ta': 'a_tibialis',
        'lf_Ta^': 'fiber_length_tibialis',
        'φ_Ta': 'pennation_angle_tibialis',
        'a_sol': 'a_soleus',
        'lf_sol^': 'fiber_length_soleus',
        'φ_sol': 'pennation_angle_soleus',
        'a_Gast': 'a_gastrocnemius',
        'lf_Gast^': 'fiber_length_gastrocnemius',
        'φ_Gast': 'pennation_angle_gastrocnemius'
    }

    # Vérifier colonnes manquantes
    missing_columns = [col for col in column_mapping.keys() if col not in df.columns]
    if missing_columns:
        raise ValueError(f"Colonnes manquantes : {missing_columns}")

    # Renommer les colonnes
    df = df.rename(columns=column_mapping)

    # Ajouter colonnes tendon_length (si pas dans Excel, initialiser à NaN)
    if 'tendon_length_tibialis' not in df.columns:
        df['tendon_length_tibialis'] = np.nan
    if 'tendon_length_soleus' not in df.columns:
        df['tendon_length_soleus'] = np.nan
    if 'tendon_length_gastrocnemius' not in df.columns:
        df['tendon_length_gastrocnemius'] = np.nan

    # Ordre des colonnes pour la sortie
    output_columns = [
        'ankle_torque',
        'q_knee', 'q_ankle',
        'a_tibialis', 'a_soleus', 'a_gastrocnemius',
        'fiber_length_tibialis', 'fiber_length_soleus', 'fiber_length_gastrocnemius',
        'pennation_angle_tibialis', 'pennation_angle_soleus', 'pennation_angle_gastrocnemius',
        'tendon_length_tibialis', 'tendon_length_soleus', 'tendon_length_gastrocnemius'
    ]

    # Extraire et transposer
    data = df[output_columns].values.T  # Shape (15, ntrials)

    # Inverser le signe du torque (première ligne) - changement de sens du repère
    data[0, :] = data[0, :]
    data[1, :] = data[1, :]
    data[2, :] = data[2, :]

    # ANGLES : deg to rad
    # q_knee (ligne 1) et q_ankle (ligne 2)
    data[1, :] = np.deg2rad(data[1, :])
    data[2, :] = np.deg2rad(data[2, :])

    # pennation angles (lignes 9, 10, 11)
    data[9, :] = np.deg2rad(data[9, :])  # pennation_angle_tibialis
    data[10, :] = np.deg2rad(data[10, :])  # pennation_angle_soleus
    data[11, :] = np.deg2rad(data[11, :])  # pennation_angle_gastrocnemius

    # LONGUEURS : cm to m (diviser par 100)
    # fiber_length (lignes 6, 7, 8)
    data[6, :] = data[6, :] / 100  # fiber_length_tibialis
    data[7, :] = data[7, :] / 100  # fiber_length_soleus
    data[8, :] = data[8, :] / 100  # fiber_length_gastrocnemius
    # tendon_length (lignes 12, 13, 14)
    data[12, :] = data[12, :] / 100  # tendon_length_tibialis
    data[13, :] = data[13, :] / 100  # tendon_length_soleus
    data[14, :] = data[14, :] / 100  # tendon_length_gastrocnemius

    print(f"✓ Shape : {data.shape}")
    print(f"  Conversions appliquées :")
    print(f"    - Torque : inversé")
    print(f"    - Angles (q, φ) : deg → rad")
    print(f"    - Longueurs (lf, lst) : cm → m")

    return data





def add_tendon_length_to_data(data, skeleton_num, casadi_function,
                              q_in_degrees=False):
    """
    Compute and fill tendon lengths in `data` using the geometric closure :
        ℓt = ℓmtu - cos(φ) * ℓm

    Vectorized over trials. Uses casadi_function['get_mtu_length'] mapped
    across all trials in a single call.

    Parameters
    ----------
    data : np.ndarray, shape (15, n_trials)
        Data matrix. Modified in place. Rows used as INPUT :
            [1]    q_knee   (rad by default)
            [2]    q_ankle  (rad by default)
            [6:9]  fiber_length (m)
            [9:12] pennation_angle (rad)
        Rows OVERWRITTEN :
            [12:15] tendon_length (m)
    skeleton_num : np.ndarray
        Musculoskeletal scalar geometric parameters (size matching the
        get_mtu_length signature).
    casadi_function : dict
        Must contain 'get_mtu_length' as a CasADi Function.
    q_in_degrees : bool, default False
        If True, q_knee and q_ankle are interpreted as degrees and converted
        to radians before computation. If False (assumes SI units), used as is.

    Returns
    -------
    data : np.ndarray
        The same array (modified in place), returned for chaining.
    """
    # Indices des lignes dans data (15, n_trials) — cohérent avec header_data
    ROW_Q_KNEE = 1
    ROW_Q_ANKLE = 2
    ROW_FIBER_LENGTH = slice(6, 9)  # tib_ant, soleus, gast
    ROW_PENNATION = slice(9, 12)
    ROW_TENDON_LENGTH = slice(12, 15)

    n_trials = data.shape[1]

    # --- Extraction et préparation des inputs --- #
    q_knee = data[ROW_Q_KNEE, :].astype(float)
    q_ankle = data[ROW_Q_ANKLE, :].astype(float)
    if q_in_degrees:
        q_knee = np.deg2rad(q_knee)
        q_ankle = np.deg2rad(q_ankle)

    fiber_length = data[ROW_FIBER_LENGTH, :].astype(float)   # (3, n_trials)
    pennation = data[ROW_PENNATION, :].astype(float)          # (3, n_trials)

    # --- Construction de l'état musculo-squelettique vectorisé --- #
    # q de taille (6, n_trials) : q[0:4]=0, q[4]=q_knee, q[5]=q_ankle
    q_full = np.zeros((6, n_trials))
    q_full[4, :] = q_knee
    q_full[5, :] = q_ankle

    # skeleton_num est constant : on le broadcast sur n_trials
    skeleton_broadcast = np.tile(np.asarray(skeleton_num).reshape(-1, 1),
                                 (1, n_trials))
    musculoskeletal_states = np.vstack([q_full, skeleton_broadcast])

    # --- Calcul vectorisé des longueurs MTU via CasADi map --- #
    # Function.map(N) crée une fonction qui applique la fonction d'origine
    # à N inputs en parallèle, en une seule évaluation.
    get_mtu_length_vec = casadi_function['get_mtu_length'].map(n_trials)
    mtu_length = np.array(get_mtu_length_vec(musculoskeletal_states))  # (3, n_trials)

    # --- Calcul des longueurs de tendon (vectorisé) --- #
    l_muscle = np.cos(pennation) * fiber_length
    tendon_length = mtu_length - l_muscle

    # --- Écriture dans data --- #
    data[ROW_TENDON_LENGTH, :] = tendon_length

    # --- Rapport console --- #
    muscle_names = ['Tibialis', 'Soleus', 'Gastrocnemius']
    print("Tendon lengths computed and added to data:")
    for i, name in enumerate(muscle_names):
        print(f"  {name:14s} : "
              "mtu length -- "
              f"min={np.nanmin(mtu_length[i]):.4f}, "
              f"max={np.nanmax(mtu_length[i]):.4f} m")
        print(f"  {name:14s} : "
              "muscle length -- "
              f"min={np.nanmin(l_muscle[i]):.4f}, "
              f"max={np.nanmax(l_muscle[i]):.4f} m")
        print(f"  {name:14s} : "
              "tendon length -- "
              f"min={np.nanmin(tendon_length[i]):.4f}, "
              f"max={np.nanmax(tendon_length[i]):.4f} m")

    return data


def get_initial_guess(muscle_tendon_parameters_num, data,
                       use_measurements=True, verbose=True):
    """
    Génère initial_guess, lower_band et upper_band avec une logique
    physiologiquement cohérente pour chaque type de paramètre.

    Stratégie :
    - lom (longueur optimale fibre) : moyenne des mesures, bornes ±25%
    - phi0 (pennation au repos) : valeur à fiber_length proche de lom, bornes ±10°
    - Fom (force max) : valeur de muscle_tendon_parameters_num, bornes ±50%
    - lst (slack length tendon) : min des mesures (tendon le moins étiré),
                                   bornes serrées ±15%

    Args:
        muscle_tendon_parameters_num : Shape (12,) - paramètres de référence
            (typiquement issus d'un modèle scalé type Rajagopal/OpenSim)
        data : Shape (15, ntrials) - données mesurées
        use_measurements : si False, utilise uniquement muscle_tendon_parameters_num
        verbose : affichage des valeurs et warnings

    Returns:
        initial_guess, upper_band, lower_band
    """

    # === Indices ===
    # data
    IDX_FL = {'ta': 6, 'sol': 7, 'gast': 8}      # fiber length
    IDX_PA = {'ta': 9, 'sol': 10, 'gast': 11}    # pennation angle
    IDX_TL = {'ta': 12, 'sol': 13, 'gast': 14}   # tendon length

    # paramètres
    IDX_LOM = {'ta': 0, 'sol': 1, 'gast': 2}
    IDX_PHI = {'ta': 3, 'sol': 4, 'gast': 5}
    IDX_FOM = {'ta': 6, 'sol': 7, 'gast': 8}
    IDX_LST = {'ta': 9, 'sol': 10, 'gast': 11}

    muscles = ['ta', 'sol', 'gast']

    # === Plages physiologiques de référence (Rajagopal 2015, Arnold 2010) ===
    # bornes anatomiques absolues (garde-fou)
    PHYSIO_BOUNDS = {
        'lom':  {'ta': (0.04, 0.10), 'sol': (0.025, 0.060), 'gast': (0.035, 0.080)},
        'phi0': {'ta': (0.05, 0.25),  'sol': (0.30, 0.65),   'gast': (0.05, 0.30)},
        'lst':  {'ta': (0.18, 0.28), 'sol': (0.20, 0.30),   'gast': (0.32, 0.45)},
    }

    initial_guess = np.zeros(12)
    lower_band = np.zeros(12)
    upper_band = np.zeros(12)

    for m in muscles:
        # --- lom : longueur optimale de fibre ---
        # Init : moyenne des fiber lengths mesurées (proxy raisonnable de lom
        # si le mouvement explore une plage autour de la longueur optimale)
        fl = data[IDX_FL[m], :]
        fl = fl[~np.isnan(fl)]
        lom_init = np.nanmean(fl) if use_measurements else muscle_tendon_parameters_num[IDX_LOM[m]]
        lom_lo = lom_init * 0.75
        lom_hi = lom_init * 1.25
        # garde-fou physiologique
        lom_lo = max(lom_lo, PHYSIO_BOUNDS['lom'][m][0])
        lom_hi = min(lom_hi, PHYSIO_BOUNDS['lom'][m][1])
        lom_init = np.clip(lom_init, lom_lo, lom_hi)

        initial_guess[IDX_LOM[m]] = lom_init
        lower_band[IDX_LOM[m]]   = lom_lo
        upper_band[IDX_LOM[m]]   = lom_hi

        # --- phi0 : pennation à la longueur optimale ---
        # On prend la pennation mesurée à fiber_length le plus proche de lom_init
        pa = data[IDX_PA[m], :]
        valid = ~(np.isnan(fl) | np.isnan(pa[:len(fl)]))
        if use_measurements and valid.any():
            idx_closest = np.argmin(np.abs(data[IDX_FL[m], :] - lom_init))
            phi_init = data[IDX_PA[m], idx_closest]
        else:
            phi_init = muscle_tendon_parameters_num[IDX_PHI[m]]
        phi_lo = phi_init - np.deg2rad(10)
        phi_hi = phi_init + np.deg2rad(10)
        phi_lo = max(phi_lo, PHYSIO_BOUNDS['phi0'][m][0])
        phi_hi = min(phi_hi, PHYSIO_BOUNDS['phi0'][m][1])
        phi_init = np.clip(phi_init, phi_lo, phi_hi)

        initial_guess[IDX_PHI[m]] = phi_init
        lower_band[IDX_PHI[m]]   = phi_lo
        upper_band[IDX_PHI[m]]   = phi_hi

        # --- Fom : force isométrique max ---
        # Pas mesurable directement, on garde la valeur scalée
        fom_init = muscle_tendon_parameters_num[IDX_FOM[m]]
        initial_guess[IDX_FOM[m]] = fom_init
        lower_band[IDX_FOM[m]]   = fom_init * 0.5
        upper_band[IDX_FOM[m]]   = fom_init * 1.5

        # --- lst : slack length tendon ---
        # Init = min des longueurs de tendon mesurées (tendon le moins étiré)
        # PAS la moyenne ! Le tendon n'est jamais plus court que lst.
        tl = data[IDX_TL[m], :]
        tl = tl[~np.isnan(tl)]
        if use_measurements and len(tl) > 0:
            lst_init = np.nanmin(tl) * 0.98  # léger margin sous le min mesuré
        else:
            lst_init = muscle_tendon_parameters_num[IDX_LST[m]]
        lst_lo = lst_init * 0.85
        lst_hi = lst_init * 1.15
        # garde-fou physiologique
        lst_lo = max(lst_lo, PHYSIO_BOUNDS['lst'][m][0])
        lst_hi = min(lst_hi, PHYSIO_BOUNDS['lst'][m][1])
        # cohérence : lst_init doit rester < min(tl) sinon contrainte tendue dès le repos
        if use_measurements and len(tl) > 0:
            lst_hi = min(lst_hi, np.nanmin(tl))
        lst_init = np.clip(lst_init, lst_lo, lst_hi)

        initial_guess[IDX_LST[m]] = lst_init
        lower_band[IDX_LST[m]]   = lst_lo
        upper_band[IDX_LST[m]]   = lst_hi

    # === Vérifications ===
    if verbose:
        print("=" * 70)
        print("Initial guess généré (avec garde-fous physiologiques)")
        print("=" * 70)
        labels = [f"{p}_{m}" for p in ['lom', 'phi0', 'Fom', 'lst'] for m in muscles]
        for i, lab in enumerate(labels):
            print(f"  {lab:12s} : init={initial_guess[i]:.4f}  "
                  f"[{lower_band[i]:.4f}, {upper_band[i]:.4f}]")

        # check sanity vs littérature
        warnings = []
        for m in muscles:
            if not (PHYSIO_BOUNDS['lom'][m][0] <= initial_guess[IDX_LOM[m]] <= PHYSIO_BOUNDS['lom'][m][1]):
                warnings.append(f"⚠ lom_{m} hors plage physiologique")
            if not (PHYSIO_BOUNDS['lst'][m][0] <= initial_guess[IDX_LST[m]] <= PHYSIO_BOUNDS['lst'][m][1]):
                warnings.append(f"⚠ lst_{m} hors plage physiologique")
        if warnings:
            print("\nAlertes :")
            for w in warnings:
                print(f"  {w}")
        print("=" * 70)

    return initial_guess, upper_band, lower_band


def plot_data(data, muscle_names=None, trial_indices=None, save_path=None):
    """
    Visualise les données expérimentales utilisées pour l'identification MTU.

    Args:
        data (np.ndarray): Shape (15, n_trials).
            Lignes :
              [0]    : torque mesuré (N.m)
              [1:3]  : angles articulaires q (rad)
              [3:6]  : activations musculaires normalisées [0, 1]
              [6:9]  : longueurs de fibres mesurées (m)
              [9:12] : angles de pennation mesurés (rad)
              [12:15]: longueurs de tendon mesurées (m)
        muscle_names (list[str], optional): noms des 3 muscles. Défaut : ['M1','M2','M3'].
        trial_indices (array-like, optional): indices x à afficher. Défaut : range(n_trials).
        save_path (str, optional): si fourni, sauvegarde la figure.

    Returns:
        fig, axes : objets matplotlib.
    """
    n_trials = data.shape[1]
    if muscle_names is None:
        muscle_names = ['M1', 'M2', 'M3']
    if trial_indices is None:
        trial_indices = np.arange(n_trials)

    # Extraction (data en shape (15, n_trials))
    torque   = data[0, :]
    q        = data[1:3, :]    # (2, n_trials)
    a        = data[3:6, :]    # (3, n_trials)
    fiber_l  = data[6:9, :]    # (3, n_trials)
    penn     = data[9:12, :]   # (3, n_trials)
    tendon_l = data[12:15, :]  # (3, n_trials)

    # Couleurs cohérentes par muscle
    colors = plt.cm.viridis(np.linspace(0.15, 0.85, 3))

    fig, axes = plt.subplots(3, 2, figsize=(13, 10), constrained_layout=True)
    fig.suptitle(f'Données expérimentales — {n_trials} trials', fontsize=13, fontweight='bold')

    # --- Torque ---
    ax = axes[0, 0]
    ax.plot(trial_indices, torque, 'o-', color='black', markersize=4, linewidth=1)
    ax.set_ylabel('Torque (N.m)')
    ax.set_title('Couple articulaire mesuré')
    ax.grid(alpha=0.3)
    ax.axhline(0, color='gray', linewidth=0.5)

    # --- Angles articulaires q ---
    ax = axes[0, 1]
    ax.plot(trial_indices, np.rad2deg(q[0, :]), 'o-', label='q1', markersize=4, linewidth=1)
    ax.plot(trial_indices, np.rad2deg(q[1, :]), 's-', label='q2', markersize=4, linewidth=1)
    ax.set_ylabel('Angle (°)')
    ax.set_title('Angles articulaires')
    ax.legend(loc='best', fontsize=9)
    ax.grid(alpha=0.3)

    # --- Activations ---
    ax = axes[1, 0]
    for m in range(3):
        ax.plot(trial_indices, a[m, :], 'o-', color=colors[m],
                label=muscle_names[m], markersize=4, linewidth=1)
    ax.set_ylabel('Activation [0–1]')
    ax.set_title('Activations musculaires')
    ax.set_ylim(-0.05, 1.05)
    ax.legend(loc='best', fontsize=9)
    ax.grid(alpha=0.3)

    # --- Longueurs de fibres ---
    ax = axes[1, 1]
    for m in range(3):
        ax.plot(trial_indices, fiber_l[m, :] * 1000, 'o-', color=colors[m],
                label=muscle_names[m], markersize=4, linewidth=1)
    ax.set_ylabel('Longueur fibre (mm)')
    ax.set_title('Longueurs de fibres mesurées')
    ax.legend(loc='best', fontsize=9)
    ax.grid(alpha=0.3)
    if (fiber_l < 0).any():
        ax.axhline(0, color='red', linewidth=0.8, linestyle='--', alpha=0.6)

    # --- Angles de pennation ---
    ax = axes[2, 0]
    for m in range(3):
        ax.plot(trial_indices, np.rad2deg(penn[m, :]), 'o-', color=colors[m],
                label=muscle_names[m], markersize=4, linewidth=1)
    ax.set_ylabel('Pennation (°)')
    ax.set_xlabel('Trial')
    ax.set_title('Angles de pennation mesurés')
    ax.legend(loc='best', fontsize=9)
    ax.grid(alpha=0.3)
    ax.axhline(0, color='gray', linewidth=0.5)

    # --- Longueurs de tendon ---
    ax = axes[2, 1]
    for m in range(3):
        ax.plot(trial_indices, tendon_l[m, :] * 1000, 'o-', color=colors[m],
                label=muscle_names[m], markersize=4, linewidth=1)
    ax.set_ylabel('Longueur tendon (mm)')
    ax.set_xlabel('Trial')
    ax.set_title('Longueurs de tendon mesurées')
    ax.legend(loc='best', fontsize=9)
    ax.grid(alpha=0.3)
    if (tendon_l < 0).any():
        ax.axhline(0, color='red', linewidth=0.8, linestyle='--', alpha=0.6)

    # --- Diagnostic console : valeurs anormales ---
    print("\n=== Diagnostic data ===")
    print(f"n_trials : {n_trials}")
    print(f"Torque         : [{torque.min():.2f}, {torque.max():.2f}] N.m")
    print(f"Fiber length   : [{fiber_l.min()*1000:.2f}, {fiber_l.max()*1000:.2f}] mm")
    print(f"Pennation      : [{np.rad2deg(penn.min()):.2f}, {np.rad2deg(penn.max()):.2f}] °")
    print(f"Tendon length  : [{tendon_l.min()*1000:.2f}, {tendon_l.max()*1000:.2f}] mm")

    neg_fl = np.argwhere(fiber_l < 0)   # (muscle, trial)
    neg_tl = np.argwhere(tendon_l < 0)
    if len(neg_fl):
        print(f"⚠ Fiber length négatif en (muscle, trial) : {neg_fl.tolist()}")
    if len(neg_tl):
        print(f"⚠ Tendon length négatif en (muscle, trial) : {neg_tl.tolist()}")
    if not len(neg_fl) and not len(neg_tl):
        print("✓ Toutes les longueurs sont positives")

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Figure sauvegardée : {save_path}")

    plt.show()
    return fig, axes

def save_muscle_tendon_parameters(muscle_tendon_parameters,
                                  osim_folder,
                                  filename='muscle_tendon_parameters_saved',
                                  muscle_names=('tibialis', 'soleus', 'gastrocnemius'),
                                  save_csv=True,
                                  save_npy=True,
                                  verbose=True):
    """
    Sauvegarde les paramètres muscle-tendon optimisés dans `osim_folder`.

    L'ordre attendu pour `muscle_tendon_parameters` (taille 12) est :
        [ℓom_TA, ℓom_SOL, ℓom_GAST,
         φo_TA,  φo_SOL,  φo_GAST,
         Fom_TA, Fom_SOL, Fom_GAST,
         ℓst_TA, ℓst_SOL, ℓst_GAST]

    Args:
        muscle_tendon_parameters (array-like, 12): paramètres optimisés.
        osim_folder (str): dossier de destination. Créé si inexistant.
        filename (str): nom de base (sans extension).
        muscle_names (tuple[str]): noms des 3 muscles dans l'ordre [TA, SOL, GAST].
        save_csv, save_npy (bool): formats à écrire.
        verbose (bool): logs.

    Returns:
        dict: chemins des fichiers écrits {'csv': ..., 'npy': ...}.
    """

    # --- Validation ---
    params = np.asarray(muscle_tendon_parameters, dtype=float).flatten()
    assert params.size == 12, f"Attendu 12 paramètres, reçu {params.size}"
    assert len(muscle_names) == 3, "3 noms de muscles attendus"

    # --- Structuration par muscle ---
    param_labels = ['optimal_fiber_length_m',
                    'pennation_angle_at_optimal_rad',
                    'max_isometric_force_N',
                    'tendon_slack_length_m']

    # Reshape : (4 paramètres, 3 muscles)
    params_matrix = params.reshape(4, 3)

    df = pd.DataFrame(params_matrix,
                      index=param_labels,
                      columns=list(muscle_names))

    # --- Préparation sortie ---
    os.makedirs(osim_folder, exist_ok=True)
    base = os.path.splitext(filename)[0]
    paths = {}

    # --- CSV : lisible pour inspection / annexes ---
    if save_csv:
        path_csv = os.path.join(osim_folder, base + '.csv')
        df_export = df.copy()
        # Ligne bonus en degrés pour lecture humaine
        df_export.loc['pennation_angle_at_optimal_deg'] = np.rad2deg(df.loc['pennation_angle_at_optimal_rad'])
        df_export.to_csv(path_csv)
        paths['csv'] = path_csv
        if verbose:
            print(f"  → CSV : {path_csv}")

    # --- NPY : recharge rapide dans le pipeline ---
    if save_npy:
        path_npy = os.path.join(osim_folder, base + '.npy')
        np.save(path_npy, params)
        paths['npy'] = path_npy
        if verbose:
            print(f"  → NPY : {path_npy}  (vecteur plat taille 12, unités SI)")

    # --- Console summary ---
    if verbose:
        print(f"\n  Timestamp : {datetime.now().isoformat(timespec='seconds')}")
        print('  Muscle-tendon parameters saved:')
        print(df.round(5).to_string())
        print()

    return paths