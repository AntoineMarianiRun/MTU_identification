"""
hot_start.py — valeurs de depart (hot start) pour un NLP d'estimation des
parametres musculo-tendineux (Hill), a partir d'essais isometriques exportes par
`save_mtu_geometry_to_xlsx` et lus avec le module de base (force_length / temp).

Deux fonctions publiques :
    hot_start_muscle(folder, name, ...)         -> {muscle: MuscleParams(f0m, lom, phi0)}
    hot_start_tendon(folder, name, f_tfl, ...)  -> {muscle: TendonParams(lst, kt, f0m)}

Muscle (isometrie, f_V = 1). Force de fibre  F_fibre = (couple/bras)/cos(pennation).
Pres de la longueur optimale, f_L ~ 1 donc F_fibre ~ F0M * a. Les essais du
"plateau" (|F_fibre|/a proche du maximum) sont ceux ou la fibre est proche de
l'optimal ; on en tire :
    f0m  = quantile(|F_fibre|/a)                 (borne haute du pic de force)
    lom  = longueur de fibre mediane du plateau  (longueur optimale)
    phi0 = pennation mediane du plateau           (angle a la longueur optimale, rad)

Tendon : on ajuste TON modele `f_tfl` (CasADi,
    f_tfl(l_tendon, f_max_tendon, t_slack_length, k_tendon) -> Force) par moindres
    carres (identify_tendon), pour (lst, kt) [+ f0m].

Separation des 3 muscles, identique pour les deux :
  * Tibial        : essais ou SEUL le tibial est actif.
  * Soleaire      : genou flechi -> gastrocnemien (bi-articulaire) detendu.
  * Gastrocnemien : genou tendu -> couple = gast + sol, on retire le soleaire.

Depend uniquement du module de base (importe ci-dessous sous force_length ;
chez toi il s'appelle peut-etre temp -> adapter la ligne d'import).
"""
import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from collections import namedtuple
import warnings


# ====================================================================== #
#  Filtres composables : fabrique -> closure (get_col, n) -> masque bool  #
#  get_col(prefix) -> (valeurs ndarray, unité str)                        #
# ====================================================================== #

def f_angle(joint, values, tol=None):
    """Ne garder que les essais dont l'angle vaut l'une des `values`.
    joint : 'knee' / 'ankle' (ou préfixe explicite). values : scalaire ou liste.
    tol   : si None -> 0.5° (ou l'équivalent en rad si fichier en SI)."""
    prefix = {'knee': 'q_knee', 'ankle': 'q_ankle'}.get(joint, joint)
    label_joint = {'q_knee': 'genou', 'q_ankle': 'cheville'}.get(prefix, prefix)
    vals = np.atleast_1d(values).astype(float)

    def inner(get_col, n):
        q, unit = get_col(prefix)
        t = tol if tol is not None else (0.5 if unit == 'deg' else np.deg2rad(0.5))
        m = np.zeros(n, dtype=bool)
        for v in vals:
            m |= np.abs(q - v) <= t
        return m

    inner.label = f"{label_joint} in {{{', '.join(f'{v:g}' for v in vals)}}}"
    return inner


def f_activation(muscle, lo=None, hi=None):
    """Garder les essais tels que lo <= activation(muscle) <= hi."""
    def inner(get_col, n):
        a, _ = get_col(f'{muscle}_activation')
        m = np.ones(n, dtype=bool)
        if lo is not None:
            m &= a >= lo
        if hi is not None:
            m &= a <= hi
        return m

    bounds = []
    if lo is not None:
        bounds.append(f'>= {lo:g}')
    if hi is not None:
        bounds.append(f'<= {hi:g}')
    inner.label = (f"a[{muscle}] " + ' & '.join(bounds)) if bounds else f"a[{muscle}] libre"
    return inner


def f_active(muscle, eps=0.05):
    """Muscle sollicité : activation >= eps."""
    f = f_activation(muscle, lo=eps)
    f.label = f"{muscle} actif (a >= {eps:g})"
    return f


def f_inactive(muscle, eps=0.05):
    """Muscle au repos : activation <= eps (exclut les essais où il est actif)."""
    f = f_activation(muscle, hi=eps)
    f.label = f"{muscle} inactif (a <= {eps:g})"
    return f


def f_where(predicate, label='critere personnalise'):
    """Échappatoire générique. predicate : `val -> masque`, où `val(prefix)`
    renvoie le tableau de valeurs d'une colonne. Permet les critères croisés :

        f_where(lambda val: (val('tibialis_activation') > 0.1)
                          & (val('soleus_activation')   < 0.05),
                label='tib actif & sol au repos')
    """
    def inner(get_col, n):
        val = lambda prefix: get_col(prefix)[0]
        return np.asarray(predicate(val), dtype=bool)

    inner.label = label
    return inner


# ====================================================================== #
#  Helpers internes partagés                                              #
# ====================================================================== #

def _load_and_filter(folder, name, filters=None, verbose=True):
    """Lit l'Excel, construit l'accès colonne, applique les filtres (ET).
    Renvoie (df, get_col, mask, n_total, labels)."""
    path = os.path.join(folder, name if name.endswith('.xlsx') else f'{name}.xlsx')
    if not os.path.exists(path):
        raise FileNotFoundError(f"Fichier Excel introuvable : {path}")
    df = pd.read_excel(path)
    n_total = len(df)

    def get_col(prefix):
        for col in df.columns:
            if col.split('[')[0] == prefix:
                m = re.search(r'\[(.*?)\]', col)
                return df[col].to_numpy(dtype=float), (m.group(1) if m else '')
        raise KeyError(f"Colonne introuvable : '{prefix}'\nColonnes : {list(df.columns)}")

    filters = list(filters) if filters else []
    sub_masks = [np.asarray(f(get_col, n_total), dtype=bool) for f in filters]
    labels = [getattr(f, 'label', 'filtre') for f in filters]

    mask = np.ones(n_total, dtype=bool)
    for sm in sub_masks:
        mask &= sm

    if verbose:
        print(f'[force_length] {mask.sum()}/{n_total} essais retenus')
        for lab, sm in zip(labels, sub_masks):
            print(f'    - {lab}: {sm.sum()}/{n_total}')

    if not mask.any():
        raise ValueError("Aucun essai ne passe tous les filtres (voir le detail par filtre).")

    return df, get_col, mask, n_total, labels


def _tendon_force(get_col, muscle, mask, torque_full):
    """force = couple / bras_de_levier_cheville [N], bras ramené en mètres."""
    ma, ma_u = get_col(f'{muscle}_moment_arm_ankle')
    ma_m = ma * (0.01 if ma_u == 'cm' else 1.0)
    ma_safe = np.where(np.abs(ma_m) > 1e-6, ma_m, np.nan)
    return (torque_full / ma_safe)[mask]


# ====================================================================== #
#  Extraction des données                                                 #
# ====================================================================== #

def extract_data(folder, name,
                 muscles=('tibialis', 'soleus', 'gastrocnemius'),
                 length='tendon_length', filters=None,
                 with_force=True, verbose=True):
    """
    Renvoie les données des essais retenus sous forme de DataFrame (mêmes
    arguments et mêmes filtres que plot_force_length, mais sans tracé).

    Colonnes :
        - partagées présentes : trial, ankle_torque, q_knee, q_ankle
        - par muscle : {muscle}_{length}, {muscle}_activation,
          {muscle}_moment_arm_ankle, et {muscle}_force[N] si with_force=True
    Les unités d'origine sont conservées dans les en-têtes.
    """
    if isinstance(muscles, str):
        muscles = [muscles]

    df, get_col, mask, n_total, _ = _load_and_filter(folder, name, filters, verbose)

    out = {}
    # colonnes partagées (uniquement si présentes dans le fichier)
    for prefix in ('trial', 'ankle_torque', 'q_knee', 'q_ankle'):
        try:
            vals, unit = get_col(prefix)
        except KeyError:
            continue
        out[prefix + (f'[{unit}]' if unit else '')] = vals[mask]

    torque, _ = get_col('ankle_torque')

    for muscle in muscles:
        ln, ln_u = get_col(f'{muscle}_{length}')
        a,  _    = get_col(f'{muscle}_activation')
        ma, ma_u = get_col(f'{muscle}_moment_arm_ankle')
        out[f'{muscle}_{length}[{ln_u}]']         = ln[mask]
        out[f'{muscle}_activation']               = a[mask]
        out[f'{muscle}_moment_arm_ankle[{ma_u}]'] = ma[mask]
        if with_force:
            out[f'{muscle}_force[N]'] = _tendon_force(get_col, muscle, mask, torque)

    return pd.DataFrame(out)


MuscleParams = namedtuple('MuscleParams', ['f0m', 'lom', 'phi0'])   # phi0 en radians
TendonParams = namedtuple('TendonParams', ['lst', 'kt', 'f0m'])


# ====================================================================== #
#  A. IDENTIFICATION TENDINEUSE (modele reel f_tfl + moindres carres)     #
# ====================================================================== #

def get_col(df, prefix):
    """Colonne d'un DataFrame commencant par `prefix` (avant '['), ramenee en SI :
    cm -> m pour longueurs/bras, le reste inchange (couple en N.m, force en N)."""
    for c in df.columns:
        if c.split('[')[0] == prefix:
            unit = c[c.find('[') + 1:c.find(']')] if '[' in c else ''
            x = df[c].to_numpy(dtype=float)
            return x * 0.01 if unit == 'cm' else x
    raise KeyError(f"Colonne '{prefix}' absente. Disponibles : {list(df.columns)}")


def _Ft(f_tfl, l, lst, kt, FoT):
    """Force tendineuse modelisee pour une longueur scalaire `l`.
    Ordre suppose : f_tfl(l_tendon, f_max_tendon, t_slack_length, k_tendon).
    Verifier avec f_tfl.name_in() ; si different, corriger ICI uniquement."""
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
    ok = np.isfinite(l) & np.isfinite(F); l, F = l[ok], F[ok]
    assert l.size > 3, f"Trop peu de points valides pour ajuster ({l.size})."
    if mode == '2p':
        FoT = float(np.max(F)) if FoT is None else FoT
        lst, kt, f0m, sol = fit_tendon_2p(f_tfl, l, F, FoT, p0, bounds)
    elif mode == '3p':
        lst, kt, f0m, sol = fit_tendon_3p(f_tfl, l, F, p0, bounds)
    else:
        raise ValueError("mode doit valoir '2p' ou '3p'.")
    return TendonParams(lst, kt, f0m), l, F, sol


def plot_tendon_fit(f_tfl, l_exp, F_exp, lst, kt, FoT, muscle_name='', ax=None):
    """Trace le modele ajuste par-dessus les donnees (+ RMSE / R2)."""
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
    ax.set_title(f'{muscle_name} - ajustement force-longueur du tendon' if muscle_name
                 else 'Ajustement force-longueur du tendon')
    ax.legend(loc='lower right', fontsize=9, framealpha=0.9); ax.grid(alpha=0.3)
    return ax, dict(rmse=rmse, r2=r2)


# ====================================================================== #
#  B. HELPERS MUSCLE (F0M, lom, phi0)                                     #
# ====================================================================== #

def _moment_arm_m(gc, muscle):
    ma, u = gc(f'{muscle}_moment_arm_ankle')
    return ma * (0.01 if u == 'cm' else 1.0)


def _pennation_rad(gc, muscle):
    """Pennation [rad] (tableau complet) ou None si absente."""
    for suffix in ('pennation_angle', 'pennation', 'alpha'):
        try:
            ang, u = gc(f'{muscle}_{suffix}')
        except KeyError:
            continue
        return np.deg2rad(ang) if u in ('deg', 'degree', 'degrees', '°') else ang
    return None


def _cos_penn(gc, muscle, mask):
    a = _pennation_rad(gc, muscle)
    if a is None:
        warnings.warn(f"[hot_start] pennation absente pour '{muscle}' -> cos = 1.")
        return np.ones(int(mask.sum()))
    return np.cos(a[mask])


def _fiber_force(gc, muscle, mask, torque):
    return _tendon_force(gc, muscle, mask, torque) / _cos_penn(gc, muscle, mask)


def _fiber_length_m(gc, muscle, mask, length):
    """Longueur de fibre [m] (cm -> m) sur le sous-ensemble masque + unite d'origine."""
    for nm in (length, 'fiber_length', 'muscle_length'):
        try:
            v, u = gc(f'{muscle}_{nm}')
            return v[mask] * (0.01 if u == 'cm' else 1.0)
        except KeyError:
            continue
    return np.full(int(mask.sum()), np.nan)


def _estimate_muscle(fiber_force, activation, fiber_length_m, pennation_rad,
                     a_min=0.3, q=0.95, plateau=0.95):
    """Estime (f0m, lom, phi0) et les masques (keep, near) des essais utilises.
      f0m  = quantile q de |F_fibre|/a  (f_L ~ 1 pres de l'optimal)
      near = essais du plateau : |F_fibre|/a >= plateau * f0m  (fibre proche de l'optimal)
      lom  = longueur de fibre mediane du plateau  [m]
      phi0 = pennation mediane du plateau           [rad]"""
    f = np.abs(np.asarray(fiber_force, float)); a = np.asarray(activation, float)
    keep = (a >= a_min) & (a > 0) & np.isfinite(f)
    if keep.sum() == 0:
        return np.nan, np.nan, np.nan, keep, keep
    ratio = np.where(keep, f / np.where(a > 0, a, np.nan), np.nan)
    f0m = float(np.quantile(ratio[keep], q))
    near = keep & (ratio >= plateau * f0m)
    if near.sum() < 3:                                   # garde-fou : top-k par ratio
        order = np.argsort(np.where(keep, ratio, -np.inf))
        k = max(3, int(np.ceil(0.05 * keep.sum())))
        near = np.zeros_like(keep); near[order[-k:]] = True
    L = np.asarray(fiber_length_m, float)
    lom = float(np.nanmedian(L[near]))
    phi0 = (float(np.nanmedian(np.asarray(pennation_rad, float)[near]))
            if pennation_rad is not None else np.nan)
    return f0m, lom, phi0, keep, near


# ====================================================================== #
#  C. HOT START MUSCLE (F0M + lom + phi0)                                 #
# ====================================================================== #

def hot_start_muscle(folder, name, *,
                     knee_ext=0.0, knee_flex=90.0,
                     eps=0.05, a_min=0.3, q=0.95, plateau=0.95,
                     length='fiber_length',
                     plot=True, save_fig=False, verbose=True):
    """Hot start de (F0M, lom, phi0) pour 'tibialis', 'soleus', 'gastrocnemius'.
    Renvoie {muscle: MuscleParams(f0m [N], lom [m], phi0 [rad])}.
    Angles genou dans les unites du fichier."""
    out, diag = {}, {}

    def _store(muscle, gc, m, Ff):
        a = gc(f'{muscle}_activation')[0][m]
        Lm = _fiber_length_m(gc, muscle, m, length)
        phi = _pennation_rad(gc, muscle); phi_m = phi[m] if phi is not None else None
        f0m, lom, phi0, keep, near = _estimate_muscle(Ff, a, Lm, phi_m, a_min, q, plateau)
        out[muscle] = MuscleParams(f0m, lom, phi0)
        diag[muscle] = dict(F=Ff, a=a, Lm=Lm, f0m=f0m, lom=lom, phi0=phi0, keep=keep, near=near)

    # 1) TIBIAL - seul actif
    filt = [f_active('tibialis', eps), f_inactive('soleus', eps), f_inactive('gastrocnemius', eps)]
    _, gc, m, _, _ = _load_and_filter(folder, name, filt, verbose)
    tau, _ = gc('ankle_torque')
    _store('tibialis', gc, m, _fiber_force(gc, 'tibialis', m, tau))

    # 2) SOLEAIRE - genou flechi (gastrocnemien detendu)
    filt = [f_angle('knee', knee_flex), f_active('soleus', eps), f_active('gastrocnemius', eps), f_inactive('tibialis', eps)]
    _, gc, m, _, _ = _load_and_filter(folder, name, filt, verbose)
    tau, _ = gc('ankle_torque')
    _store('soleus', gc, m, _fiber_force(gc, 'soleus', m, tau))

    # 3) GASTROCNEMIEN - genou tendu ; on retire le soleaire (via son F0M)
    filt = [f_angle('knee', knee_ext), f_active('soleus', eps), f_active('gastrocnemius', eps), f_inactive('tibialis', eps)]
    _, gc, m, _, _ = _load_and_filter(folder, name, filt, verbose)
    tau, _ = gc('ankle_torque'); a_pf = gc('soleus_activation')[0]
    r_sol = _moment_arm_m(gc, 'soleus')[m]; r_gas = _moment_arm_m(gc, 'gastrocnemius')[m]
    tau_sol = out['soleus'].f0m * a_pf[m] * _cos_penn(gc, 'soleus', m) * r_sol
    r_gas_safe = np.where(np.abs(r_gas) > 1e-6, r_gas, np.nan)
    Ff_gas = ((tau[m] - tau_sol) / r_gas_safe) / _cos_penn(gc, 'gastrocnemius', m)
    _store('gastrocnemius', gc, m, Ff_gas)

    if verbose:
        print('\n[hot_start] muscle (F0M | lom | phi0) :')
        for k, p in out.items():
            print(f'    {k:14s}: F0M = {p.f0m:7.1f} N | lom = {p.lom*100:5.2f} cm | phi0 = {np.rad2deg(p.phi0):4.1f} deg')
    if plot:
        _plot_muscle(diag, out, folder, name, save_fig)
    return out


def _plot_muscle(diag, out, folder, name, save_fig):
    muscles = ['tibialis', 'soleus', 'gastrocnemius']
    titles = {'tibialis': 'Tibial anterieur (seul actif)', 'soleus': 'Soleaire (genou flechi)',
              'gastrocnemius': 'Gastrocnemien (genou tendu - soleaire)'}
    fig, axes = plt.subplots(2, 2, figsize=(11, 8.4), constrained_layout=True); sc = None
    for ax, mus in zip(axes.flat[:3], muscles):
        d = diag[mus]; x, f, a = d['Lm'] * 100, np.abs(d['F']), d['a']
        keep, near = d['keep'], d['near']; p = out[mus]
        ax.scatter(x[~keep], f[~keep], c='0.85', s=14, edgecolors='none', zorder=1, label='a < a_min')
        sc = ax.scatter(x[keep], f[keep], c=a[keep], cmap='RdYlBu_r', vmin=0, vmax=1,
                        s=22 + a[keep] * 70, alpha=0.85, edgecolors='none', zorder=2)
        ax.scatter(x[near], f[near], s=90, facecolors='none', edgecolors='k', linewidths=1.1,
                   zorder=3, label='plateau (lom, phi0)')
        ax.axhline(p.f0m, color='k', ls='--', lw=1.2, label=f"F0M ~ {p.f0m:.0f} N")
        ax.axvline(p.lom * 100, color='0.35', ls=':', lw=1.4, label=f"lom ~ {p.lom*100:.2f} cm")
        ax.text(0.03, 0.05, f"phi0 ~ {np.rad2deg(p.phi0):.1f} deg", transform=ax.transAxes,
                fontsize=9, bbox=dict(boxstyle='round', fc='white', ec='0.6', alpha=0.9))
        ax.set_title(titles[mus], fontsize=10); ax.set_xlabel('longueur de fibre (cm)')
        ax.set_ylabel('Force de fibre (N)'); ax.grid(True, lw=0.4, alpha=0.4); ax.legend(fontsize=7.5, loc='upper left')
    if sc is not None:
        fig.colorbar(sc, ax=list(axes.flat[:3]), label='Activation', shrink=0.85)
    axt = axes.flat[3]; axt.axis('off')
    header = ['muscle', 'F0M [N]', 'lom [cm]', 'phi0 [deg]']
    rows = [[mus, f"{out[mus].f0m:.0f}", f"{out[mus].lom*100:.2f}", f"{np.rad2deg(out[mus].phi0):.1f}"]
            for mus in muscles]
    tbl = axt.table(cellText=rows, colLabels=header, loc='center', cellLoc='center')
    tbl.auto_set_font_size(False); tbl.set_fontsize(9); tbl.scale(1, 1.7)
    axt.set_title('Hot start muscle', fontsize=10)
    fig.suptitle('Hot start muscle - F0M, lom (longueur optimale), phi0 (pennation a l\'optimal)', fontsize=11)
    if save_fig:
        import os
        out_png = os.path.join(folder or '.', f'{os.path.splitext(name)[0]}_hot_start_muscle.png')
        fig.savefig(out_png, dpi=150, bbox_inches='tight'); print(f'[hot_start] figure : {out_png}')
    plt.show()


# ====================================================================== #
#  D. HOT START TENDON (utilise ton f_tfl + identify_tendon)              #
# ====================================================================== #

def hot_start_tendon(folder, name, f_tfl, *,
                     ankle=(-10, 0), knee_flex=90, knee_ext=0,
                     mode='2p', f0m=None, plot=True, verbose=False):
    """Hot start des parametres tendineux (lst, kt) [+ f0m] pour soleaire, tibial
    anterieur et gastrocnemien, en ajustant TON modele `f_tfl` par moindres carres
    (identify_tendon), restreint aux angles de cheville `ankle` (defaut -10 et 0).

    `f0m` : dict {muscle: F0M} OU {muscle: MuscleParams} (ex. sortie de
            hot_start_muscle) OU None -> estime en interne (sans restriction cheville).
    mode='2p' : f0m fige -> ajuste (lst, kt) ; mode='3p' : (lst, kt, f0m) libres.
    Renvoie {muscle: TendonParams(lst, kt, f0m)}."""
    if f0m is None:
        mp = hot_start_muscle(folder, name, knee_ext=knee_ext, knee_flex=knee_flex,
                              plot=False, verbose=False)
        f0m = {k: v.f0m for k, v in mp.items()}
    else:
        f0m = {k: (v.f0m if hasattr(v, 'f0m') else float(v)) for k, v in f0m.items()}

    qa = [f_angle('ankle', list(ankle))] if ankle is not None else []
    out, packed = {}, {}

    # 1) SOLEAIRE - genou flechi (gastrocnemien detendu)
    d = extract_data(folder, name, muscles=['soleus'], length='tendon_length',
                     filters=qa + [f_angle('knee', knee_flex), f_inactive('tibialis')], verbose=verbose)
    p_sol, l_sol, F_sol, _ = identify_tendon(
        f_tfl, get_col(d, 'soleus_tendon_length'), get_col(d, 'soleus_force'),
        mode=mode, FoT=f0m.get('soleus'))
    out['soleus'] = p_sol; packed['soleus'] = (l_sol, F_sol)

    # 2) TIBIAL ANTERIEUR - seul dorsiflechisseur actif (|.| : couple oppose)
    d = extract_data(folder, name, muscles=['tibialis'], length='tendon_length',
                     filters=qa + [ f_inactive('soleus'), f_inactive('gastrocnemius')], verbose=verbose)
    p_tib, l_tib, F_tib, _ = identify_tendon(
        f_tfl, get_col(d, 'tibialis_tendon_length'), np.abs(get_col(d, 'tibialis_force')),
        mode=mode, FoT=f0m.get('tibialis'))
    out['tibialis'] = p_tib; packed['tibialis'] = (l_tib, F_tib)

    # 3) GASTROCNEMIEN - genou tendu ; soleaire retire via SA courbe f_tfl (etape 1)
    d = extract_data(folder, name, muscles=['soleus', 'gastrocnemius'], length='tendon_length',
                     filters=qa + [f_angle('knee', knee_ext), f_inactive('tibialis')], verbose=verbose)

    tau_total = get_col(d, 'ankle_torque')
    l_sol_g   = get_col(d, 'soleus_tendon_length')
    ma_sol_g  = get_col(d, 'soleus_moment_arm_ankle')
    l_gas_raw = get_col(d, 'gastrocnemius_tendon_length')
    ma_gas    = get_col(d, 'gastrocnemius_moment_arm_ankle')
    F_sol_est = np.array([_Ft(f_tfl, l, p_sol.lst, p_sol.kt, p_sol.f0m) for l in l_sol_g])
    F_gas_raw = (tau_total - F_sol_est * ma_sol_g) / ma_gas
    valid = np.isfinite(F_gas_raw) & (F_gas_raw > 0)
    if valid.sum() < F_gas_raw.size:
        print(f"[Gastrocnem.] {F_gas_raw.size - valid.sum()}/{F_gas_raw.size} essais ecartes (F_gastro <= 0).")
    p_gas, l_gas, F_gas, _ = identify_tendon(
        f_tfl, l_gas_raw[valid], F_gas_raw[valid], mode=mode, FoT=f0m.get('gastrocnemius'))
    out['gastrocnemius'] = p_gas; packed['gastrocnemius'] = (l_gas, F_gas)

    print("\n=== Hot start tendon (cheville =", ankle, "| mode =", mode, ") ===")
    for nm, key in [('Soleaire', 'soleus'), ('Tibial ant.', 'tibialis'), ('Gastrocnem.', 'gastrocnemius')]:
        p = out[key]
        print(f"{nm:13s}: lst = {p.lst*100:6.2f} cm | kt = {p.kt:6.2f} | f0m = {p.f0m:8.0f} N")

    if plot:
        fig, axes = plt.subplots(1, 3, figsize=(18, 5))
        for ax, key, nm in zip(axes, ['soleus', 'tibialis', 'gastrocnemius'],
                               ['Soleaire', 'Tibial anterieur', 'Gastrocnemien']):
            l, F = packed[key]; p = out[key]
            plot_tendon_fit(f_tfl, l, F, p.lst, p.kt, p.f0m, muscle_name=nm, ax=ax)
        fig.tight_layout(); plt.show()
    return out


def merge_tendon_muscle_param(mus,ten):
    muscles = ['tibialis', 'soleus', 'gastrocnemius']  # ton ordre de muscles
    order = ['l0m', 'phi0', 'f0m', 'lst', 'kt']  # ton ordre de parametres

    def value(param, m):
        return {
            'l0m': mus[m].lom,
            'phi0': mus[m].phi0,
            'f0m': mus[m].f0m,
            'lst': ten[m].lst,
            'kt': ten[m].kt,
        }[param]

    # regroupe par PARAMETRE : [l0m(tib,sol,gas), phi0(tib,sol,gas), f0m(...), lst(...), km(...), kt(...)]
    x0 = np.concatenate([[value(p, m) for m in muscles] for p in order])

    return x0