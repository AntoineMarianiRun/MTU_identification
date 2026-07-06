# cette file .py est la pour determiner les parametres du tendon
import manipfun
import import_functions
import useful
import temp
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from collections import namedtuple

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
osim_path, train_path, test_path = manipfun.select_folder_and_get_files()  # get the interest folder

# 6.1.3 name and path of the files
osim_folder, osim_name = manipfun.split_path_name(osim_path)
train_folder, train_name = manipfun.split_path_name(train_path)
test_folder, test_name = manipfun.split_path_name(test_path)

# 6.1.4 import the osim scaled as a python variable
mtu_params = ['l0m', 'phi0', 'f0m', 'lst', 'kt']
param_config = {
    'l0m': 'sym',
    'phi0': 'sym',
    'f0m': 'sym',
    'km': 'fixed',
    'lst': 'sym',
    'kt': 'sym',
}

skeleton_num, muscle_tendon_parameters_num = import_functions.get_model_osim_scaled(
    osim_folder, osim_name, mtu_params=mtu_params)
casadi_function, unknown_parameters, definition = useful.get_model_equation(param_config=param_config)

# 6.1.5 training data set and import test data set
data_train = useful.load_data_from_xlsx(train_folder, train_name, header)
data_test = useful.load_data_from_xlsx(test_folder, test_name, header)

# 6.1.6 add tendon length to data set and import test data set
data_train = manipfun.add_tendon_length_to_data(data_train, skeleton_num, casadi_function)
data_test = manipfun.add_tendon_length_to_data(data_test, skeleton_num, casadi_function)

useful.save_mtu_geometry_to_xlsx(data_train, skeleton_num, casadi_function,
                                 train_folder, 'geometrie_train_v2',
                                 convert_units=True, q_in_degrees=False)
file_name = 'geometrie_train_v2'

# fonction tendon commune aux 3 muscles :
#   f(l_tendon, f_max_tendon, t_slack_length, k_tendon) -> Force_tendon
f_tfl = casadi_function['tendon_force_length_single_muscle']


# ====================================================================== #
#  OUTILS GÉNÉRIQUES D'IDENTIFICATION DES PARAMÈTRES TENDINEUX            #
# ====================================================================== #
TendonParams = namedtuple('TendonParams', ['lst', 'kt', 'f0m'])


def get_col(df, prefix):
    """Colonne du DataFrame commençant par `prefix` (avant le '['), ramenée en SI :
    cm -> m pour les longueurs et bras de levier, le reste inchange (couple en N.m)."""
    for c in df.columns:
        if c.split('[')[0] == prefix:
            unit = c[c.find('[') + 1:c.find(']')] if '[' in c else ''
            x = df[c].to_numpy(dtype=float)
            return x * 0.01 if unit == 'cm' else x
    raise KeyError(f"Colonne '{prefix}' absente. Disponibles : {list(df.columns)}")


def _Ft(f_tfl, l, lst, kt, FoT):
    """Force tendineuse modelisee pour une longueur scalaire `l`.
    /!\\ Ordre des entrees suppose : f_tfl(l_tendon, f_max_tendon, t_slack_length, k_tendon).
        Verifier avec f_tfl.name_in() ; si different, corriger l'appel ICI uniquement."""
    return float(f_tfl(l, FoT, lst, kt))


def fit_tendon_2p(f_tfl, l_exp, F_exp, FoT, p0=None, bounds=None):
    """Ajuste (lst, kt), f0m fige a FoT. Renvoie (lst, kt, FoT, sol)."""
    lst_min = float(np.min(l_exp))
    if p0 is None:     p0 = [0.95 * lst_min, 35.0]
    if bounds is None: bounds = ([0.5 * lst_min, 10.0], [1.5 * lst_min, 100.0])

    def res(p):
        lst, kt = p
        return np.array([_Ft(f_tfl, l, lst, kt, FoT) for l in l_exp]) - F_exp

    sol = least_squares(res, x0=p0, bounds=bounds, method='trf',
                        x_scale='jac', ftol=1e-12, xtol=1e-12)
    lst, kt = sol.x
    return lst, kt, FoT, sol


def fit_tendon_3p(f_tfl, l_exp, F_exp, p0=None, bounds=None):
    """Ajuste (lst, kt, f0m) tous libres. Renvoie (lst, kt, f0m, sol)."""
    lst_min = float(np.min(l_exp)); F0 = float(np.max(F_exp))
    if p0 is None:     p0 = [0.95 * lst_min, 35.0, F0]
    if bounds is None: bounds = ([0.5 * lst_min, 10.0, 0.1 * F0],
                                 [1.5 * lst_min, 100.0, 5.0 * F0])

    def res(p):
        lst, kt, FoT = p
        return np.array([_Ft(f_tfl, l, lst, kt, FoT) for l in l_exp]) - F_exp

    sol = least_squares(res, x0=p0, bounds=bounds, method='trf',
                        x_scale='jac', ftol=1e-12, xtol=1e-12)
    lst, kt, FoT = sol.x
    return lst, kt, FoT, sol


def identify_tendon(f_tfl, l_exp, F_exp, mode='3p', FoT=None, p0=None, bounds=None):
    """Nettoie les NaN, ajuste le modele tendineux et renvoie
    (TendonParams(lst, kt, f0m), l_clean, F_clean, sol).
    mode='2p' : f0m fige (a FoT, sinon max(F_exp)) ; mode='3p' : f0m libre."""
    l = np.asarray(l_exp, float); F = np.asarray(F_exp, float)
    ok = np.isfinite(l) & np.isfinite(F)
    l, F = l[ok], F[ok]
    assert l.size > 3, f"Trop peu de points valides pour ajuster les parametres ({l.size})."
    if mode == '2p':
        FoT = float(np.max(F)) if FoT is None else FoT
        lst, kt, f0m, sol = fit_tendon_2p(f_tfl, l, F, FoT, p0, bounds)
    elif mode == '3p':
        lst, kt, f0m, sol = fit_tendon_3p(f_tfl, l, F, p0, bounds)
    else:
        raise ValueError("mode doit valoir '2p' ou '3p'.")
    return TendonParams(lst, kt, f0m), l, F, sol


def plot_tendon_fit(f_tfl, l_exp, F_exp, lst, kt, FoT, muscle_name='', ax=None):
    """Trace le modele ajuste par-dessus les donnees experimentales (+ RMSE / R2)."""
    l_grid  = np.linspace(min(l_exp.min(), lst), l_exp.max() * 1.01, 300)
    F_curve = np.array([_Ft(f_tfl, l, lst, kt, FoT) for l in l_grid])
    F_pred  = np.array([_Ft(f_tfl, l, lst, kt, FoT) for l in l_exp])

    rmse = float(np.sqrt(np.mean((F_pred - F_exp) ** 2)))
    ss_res, ss_tot = np.sum((F_exp - F_pred) ** 2), np.sum((F_exp - F_exp.mean()) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float('nan')

    if ax is None:
        _, ax = plt.subplots(figsize=(7, 5))
    ax.scatter(l_exp * 100, F_exp, s=28, color='#1f77b4', alpha=0.75,
               edgecolor='white', linewidth=0.5, label='Donnees experimentales', zorder=3)
    ax.plot(l_grid * 100, F_curve, color='#d62728', lw=2.2, label='Modele ajuste', zorder=2)
    ax.axvline(lst * 100, ls='--', color='gray', lw=1, label=fr'$l_{{st}}$ = {lst*100:.2f} cm')
    ax.text(0.03, 0.97,
            fr'$l_{{st}}$ = {lst*100:.2f} cm''\n'fr'$k_t$ = {kt:.1f}''\n'
            fr'$F_0^M$ = {FoT:.0f} N''\n'fr'RMSE = {rmse:.1f} N   $R^2$ = {r2:.4f}',
            transform=ax.transAxes, va='top', fontsize=9,
            bbox=dict(boxstyle='round', fc='white', ec='0.6', alpha=0.9))
    ax.set_xlabel(r'Longueur du tendon $l_t$ [cm]')
    ax.set_ylabel(r'Force tendineuse $F_t$ [N]')
    title = f'{muscle_name} - ajustement force-longueur du tendon' if muscle_name \
            else 'Ajustement force-longueur du tendon'
    ax.set_title(title)
    ax.legend(loc='lower right', fontsize=9, framealpha=0.9); ax.grid(alpha=0.3)
    return ax, dict(rmse=rmse, r2=r2)


MODE = '2p'   # '3p' : (lst, kt, f0m) libres ; '2p' : f0m fige a la force max mesuree

########################################################################################################################
#  1. SOLEAIRE - genou flechi (90 deg) : seul plantarflechisseur actif
#     (gastrocnemien detendu, tibial inactif)
########################################################################################################################
data_soleus = temp.extract_data(train_folder, file_name,
    muscles=['soleus'], length='tendon_length',
    filters=[
        temp.f_angle('knee', 90),
        temp.f_angle('ankle', [0]),
        temp.f_inactive('tibialis'),
    ])

p_sol, l_sol, F_sol, sol_sol = identify_tendon(
    f_tfl,
    get_col(data_soleus, 'soleus_tendon_length'),   # m
    get_col(data_soleus, 'soleus_force'),           # N  (= couple / bras soleaire)
    mode=MODE)
print(f"[Soleaire]    succes={sol_sol.success} | lst={p_sol.lst*100:.2f} cm | "
      f"kt={p_sol.kt:.2f} | f0m={p_sol.f0m:.0f} N")

########################################################################################################################
#  2. TIBIAL ANTERIEUR - meme methode : seul dorsiflechisseur actif
#     (soleaire ET gastrocnemien inactifs)
########################################################################################################################
#  /!\ Adapter les filtres a ton protocole (angles, seuils d'activation).
#  /!\ Convention de signe : si la force tibiale ressort negative (couple de dorsiflexion
#      oppose aux plantaires), il faut soit l'option absolute de extract_data, soit np.abs.
data_tibialis = temp.extract_data(train_folder, file_name,
    muscles=['tibialis'], length='tendon_length',
    filters=[
        temp.f_angle('knee', 0),
        temp.f_angle('ankle', [-10]),
        temp.f_inactive('soleus'),
        temp.f_inactive('gastrocnemius'),
    ])

p_tib, l_tib, F_tib, sol_tib = identify_tendon(
    f_tfl,
    get_col(data_tibialis, 'tibialis_tendon_length'),  # m
    get_col(data_tibialis, 'tibialis_force'),          # N  (= couple / bras tibial)
    mode=MODE)
print(f"[Tibial ant.] succes={sol_tib.success} | lst={p_tib.lst*100:.2f} cm | "
      f"kt={p_tib.kt:.2f} | f0m={p_tib.f0m:.0f} N")

########################################################################################################################
#  3. GASTROCNEMIEN - genou tendu (0 deg) : soleaire + gastrocnemien actifs
#     On retire la part soleaire estimee via SA courbe force-longueur :
#         couple_gastro = couple_total - F_soleaire(l_t_soleaire) * bras_soleaire
#         F_gastro      = couple_gastro / bras_gastro
########################################################################################################################
data_gastro = temp.extract_data(train_folder, file_name,
    muscles=['soleus', 'gastrocnemius'], length='tendon_length',
    filters=[
        temp.f_angle('knee', [0]),
        temp.f_angle('ankle', [-10]),
        temp.f_inactive('tibialis'),
    ])

# colonnes brutes (memes essais, lignes alignees)
tau_total = get_col(data_gastro, 'ankle_torque')                    # N.m  (couple total cheville)
l_sol_g   = get_col(data_gastro, 'soleus_tendon_length')            # m
ma_sol_g  = get_col(data_gastro, 'soleus_moment_arm_ankle')         # m
l_gas_raw = get_col(data_gastro, 'gastrocnemius_tendon_length')     # m
ma_gas    = get_col(data_gastro, 'gastrocnemius_moment_arm_ankle')  # m

# part soleaire : force tendineuse estimee depuis SA longueur mesuree (parametres de l'etape 1)
F_sol_est = np.array([_Ft(f_tfl, l, p_sol.lst, p_sol.kt, p_sol.f0m) for l in l_sol_g])  # N
tau_sol   = F_sol_est * ma_sol_g                                   # N.m
tau_gas   = tau_total - tau_sol                                    # N.m  (part gastrocnemienne)
F_gas_raw = tau_gas / ma_gas                                       # N

# garde les essais physiquement plausibles (force gastrocnemienne positive)
valid = np.isfinite(F_gas_raw) & (F_gas_raw > 0)
if valid.sum() < F_gas_raw.size:
    print(f"[Gastrocnem.] {F_gas_raw.size - valid.sum()}/{F_gas_raw.size} essais ecartes "
          f"(F_gastro <= 0 -> verifier conventions de signe couple/bras).")

p_gas, l_gas, F_gas, sol_gas = identify_tendon(
    f_tfl, l_gas_raw[valid], F_gas_raw[valid], mode=MODE)
print(f"[Gastrocnem.] succes={sol_gas.success} | lst={p_gas.lst*100:.2f} cm | "
      f"kt={p_gas.kt:.2f} | f0m={p_gas.f0m:.0f} N")

########################################################################################################################
#  Figure recapitulative (3 panneaux) + tableau de synthese
########################################################################################################################
fig, axes = plt.subplots(1, 3, figsize=(18, 5))
plot_tendon_fit(f_tfl, l_sol, F_sol, *p_sol, muscle_name='Soleaire',         ax=axes[0])
plot_tendon_fit(f_tfl, l_tib, F_tib, *p_tib, muscle_name='Tibial anterieur', ax=axes[1])
plot_tendon_fit(f_tfl, l_gas, F_gas, *p_gas, muscle_name='Gastrocnemien',    ax=axes[2])
fig.tight_layout()
plt.show()

print("\n=== Parametres tendineux identifies ===")
for name, p, s in [('Soleaire', p_sol, sol_sol),
                   ('Tibial ant.', p_tib, sol_tib),
                   ('Gastrocnem.', p_gas, sol_gas)]:
    print(f"{name:13s}: lst = {p.lst*100:6.2f} cm | kt = {p.kt:6.2f} | "
          f"f0m = {p.f0m:8.0f} N | succes = {s.success}")