import numpy as np
from useful import solve_all_muscles
import matplotlib
matplotlib.use('TkAgg')  # backend interactif
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.widgets import Slider, RadioButtons, Button


def plotmodel(ax, muscle_origins, muscle_insertions, joint_centers, via_point):
    ax.clear()
    ax.set_title("Model", fontsize=14)

    origins = np.array(muscle_origins)
    insertions = np.array(muscle_insertions)
    markers = np.array(joint_centers)
    via_point = np.array(via_point)

    ax.scatter(origins[0], origins[1], origins[2], c='r', label='muscle origin', marker='o')
    ax.scatter(insertions[0], insertions[1], insertions[2], c='b', label='muscle insertion', marker='o')
    ax.scatter(markers[0], markers[1], markers[2], c='k', label='joint center', marker='o')
    ax.scatter(via_point[0], via_point[1], via_point[2], c='y', label='tibialis via point', marker='o')

    for i in range(markers.shape[1] - 1):
        ax.plot([markers[0, i], markers[0, i + 1]],
                [markers[1, i], markers[1, i + 1]],
                [markers[2, i], markers[2, i + 1]], color='black')

    for i in range(1, origins.shape[1]):
        ax.plot([origins[0, i], insertions[0, i]],
                [origins[1, i], insertions[1, i]],
                [origins[2, i], insertions[2, i]], color='red')

    ax.plot([origins[0, 0], via_point[0, 0]],
            [origins[1, 0], via_point[1, 0]],
            [origins[2, 0], via_point[2, 0]], color='red')

    ax.plot([via_point[0, 0], insertions[0, 0]],
            [via_point[1, 0], insertions[1, 0]],
            [via_point[2, 0], insertions[2, 0]], color='red')

    ax.legend()
    ax.set_xlim([-1, 1])
    ax.set_ylim([-1, 1])
    ax.set_zlim([-1, 1])
    ax.set_box_aspect([1, 1, 1])
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

def interactive_model(skeleton_num, muscle_tendon_parameters_num, casadi_function):
    """
    Plot interactif du modèle musculo-squelettique.

    Nouveautés par rapport à la version précédente :
      - les 18 paramètres muscle-tendon [l0m, phi0, f0m, km, lst, kt] x 3 muscles
        sont modifiables via des sliders REPLIABLES PAR MUSCLE (un RadioButtons
        sélectionne le muscle actif ; 6 sliders affichés à la fois).
      - tous les anciens `print` sont supprimés : les résultats sont écrits dans
        une SECONDE FIGURE, en texte monospace sélectionnable / copiable
        (clic-glisser sur le texte puis Ctrl+C avec le backend Qt/Tk).

    Convention de rangement de `muscle_tendon_parameters_num` : 18x1 APLATI PAR
    PARAMETRE -> [l0m_t, l0m_s, l0m_g, phi0_t, phi0_s, phi0_g, f0m_t, ...].
    Le paramètre p (0..5) du muscle m (0..2) est donc à l'index  p*3 + m.
    """

    n_musc = 3
    n_par = 6
    muscle_names = ['tibialis', 'soleus', 'gastrocnemius']
    short_names = ['tibialis', 'soleus', 'gastroc']
    param_names = ['l0m', 'phi0', 'f0m', 'km', 'lst', 'kt']
    # paramètres exprimés en degrés à l'affichage (phi0)
    param_is_angle = [False, True, False, False, False, False]

    q_init = np.zeros(6)
    muscle_tendon_parameters_num = np.array(muscle_tendon_parameters_num, dtype=float).flatten()
    assert muscle_tendon_parameters_num.size == n_musc * n_par, \
        f'attendu {n_musc * n_par} paramètres, reçu {muscle_tendon_parameters_num.size}'

    # copie de travail mutée par les sliders (toujours rangée par paramètre)
    mtu_params = muscle_tendon_parameters_num.copy()

    def pidx(p, m):
        """index aplati (rangé par paramètre) du paramètre p du muscle m."""
        return p * n_musc + m

    # =====================================================================
    # FIGURE « Modèle interactif » : vue 3D (moitié haute)
    #                                + architecture musculaire (moitié basse)
    # =====================================================================
    fig = plt.figure(figsize=(11, 11), num='Modèle interactif')

    # moitié supérieure : vue 3D
    gs_top = GridSpec(
        1, 1, figure=fig,
        left=0.05, right=0.95, bottom=0.52, top=0.97,
    )
    ax = fig.add_subplot(gs_top[0, 0], projection='3d')

    arch_colors = ['#d62728', '#2ca02c', '#9467bd']

    # moitié inférieure : 3 diagrammes d'architecture côte à côte
    gs_bot = GridSpec(
        1, 3, figure=fig,
        left=0.04, right=0.97, bottom=0.04, top=0.46, wspace=0.15,
    )
    ax_arch = [fig.add_subplot(gs_bot[0, i]) for i in range(3)]

    # =====================================================================
    # FIGURE « Contrôles » : tous les sliders (zones A / B / C)
    # =====================================================================
    fig_ctrl = plt.figure(figsize=(7, 9), num='Contrôles')

    # =====================================================================
    # FIGURE « Résultats » : panneau texte + bouton « Copier »
    # =====================================================================
    fig_txt = plt.figure(figsize=(7.5, 9), num='Résultats')
    ax_txt = fig_txt.add_axes([0.0, 0.06, 1.0, 0.93])
    ax_txt.axis('off')
    txt_handle = ax_txt.text(
        0.02, 0.98, '', transform=ax_txt.transAxes,
        va='top', ha='left', family='monospace', fontsize=8.5,
    )

    # stocke le dernier rapport en clair pour la copie
    report_store = {'text': ''}

    def copy_to_clipboard(event=None):
        text = report_store['text']
        # 1) tentative via le clipboard Tk du backend matplotlib
        try:
            tkwidget = fig_txt.canvas.get_tk_widget()
            tkwidget.clipboard_clear()
            tkwidget.clipboard_append(text)
            tkwidget.update()  # garde le contenu après fermeture de l'event loop
            btn_copy.label.set_text('Copié ✓')
            fig_txt.canvas.draw_idle()
            return
        except Exception:
            pass
        # 2) repli pyperclip (Qt / autres backends)
        try:
            import pyperclip
            pyperclip.copy(text)
            btn_copy.label.set_text('Copié ✓')
            fig_txt.canvas.draw_idle()
        except Exception:
            btn_copy.label.set_text('Copie indispo')
            fig_txt.canvas.draw_idle()

    ax_btn = fig_txt.add_axes([0.35, 0.01, 0.30, 0.04])
    btn_copy = Button(ax_btn, 'Copier les résultats')
    btn_copy.on_clicked(copy_to_clipboard)

    def set_report(text):
        report_store['text'] = text
        txt_handle.set_text(text)
        btn_copy.label.set_text('Copier les résultats')  # reset du libellé
        fig_txt.canvas.draw_idle()

    # ---------------------------------------------------------------
    # helper : schéma vue de dessus d'un muscle penné
    # ---------------------------------------------------------------
    def draw_muscle_architecture(axm, l_f, phi, l_t, name, color):
        axm.clear()
        l_f = float(np.array(l_f).flatten()[0])
        phi = float(np.array(phi).flatten()[0])
        l_t = float(np.array(l_t).flatten()[0])

        h = l_f * np.sin(phi)
        dx = l_f * np.cos(phi)
        n_fasc = 5
        apo_length = n_fasc * dx

        axm.plot([0, apo_length], [0, 0], 'k-', lw=2)
        axm.plot([dx, apo_length + dx], [h, h], 'k-', lw=2)
        for i in range(n_fasc + 1):
            x0 = i * dx
            axm.plot([x0, x0 + dx], [0, h], color=color, lw=1.3, alpha=0.85)
        axm.plot([apo_length + dx, apo_length + dx + l_t], [h, h], color='#1f77b4', lw=3.5)

        r = 0.35 * dx
        theta = np.linspace(0, phi, 30)
        axm.plot(r * np.cos(theta), r * np.sin(theta), 'k-', lw=0.8)
        axm.text(r * 1.25 * np.cos(phi / 2), r * 1.25 * np.sin(phi / 2), r'$\varphi$', fontsize=8)

        axm.set_title(
            f'{name}  |  $l_f$={l_f * 100:.2f} cm  |  '
            f'$\\varphi$={np.rad2deg(phi):.1f}°  |  $l_t$={l_t * 100:.2f} cm',
            fontsize=9,
        )
        axm.set_aspect('equal')
        axm.axis('off')

    # ---------------------------------------------------------------
    # helpers d'extraction "best effort" depuis results[muscle]
    # ---------------------------------------------------------------
    def fmt_residual(res):
        # res['residuals'] : DM vecteur (CasADi) -> ||·|| + composantes
        val = res.get('residuals', None)
        if val is None:
            return 'n/a'
        arr = np.array(val).flatten().astype(float)
        comp = ', '.join(f'{v:.2e}' for v in arr)
        return f'||{np.linalg.norm(arr):.2e}||  [{comp}]'

    def fmt_success(res):
        # res['state'] : chaîne (ex. 'success')
        return str(res.get('state', 'n/a'))

    # ---------------------------------------------------------------
    # callback principal
    # ---------------------------------------------------------------
    def update_plot(val=None):
        a = np.array([muscle_sliders[i].val for i in range(3)])
        q = np.array([slider_q[i].val for i in range(6)])
        q[3:] = np.deg2rad(q[3:])

        musculoskeletal_states_num = np.concatenate((q, skeleton_num))
        neuromusculoskeletal_state_num = np.concatenate([a, musculoskeletal_states_num])

        origins, insertions, via_point, markers = casadi_function['forward_kinematics'](musculoskeletal_states_num)
        plotmodel(ax, origins, insertions, markers, via_point)

        mtu_length = casadi_function['get_mtu_length'](musculoskeletal_states_num)
        mtu_moment_arm = casadi_function['get_moment_arm'](musculoskeletal_states_num)

        # === équilibre musculaire (avec les paramètres COURANTS) ===
        results = solve_all_muscles(a, mtu_length, mtu_params, casadi_function)

        x_opt = {
            'tibialis': results['tibialis']['x_opt'],
            'soleus': results['soleus']['x_opt'],
            'gastrocnemius': results['gastrocnemius']['x_opt'],
        }
        # rooted_variables = [fiber length, pennation angle, tendon length]
        rooted_variables = np.array([
            x_opt['tibialis'][0, 0], x_opt['soleus'][0, 0], x_opt['gastrocnemius'][0, 0],
            x_opt['tibialis'][1, 0], x_opt['soleus'][1, 0], x_opt['gastrocnemius'][1, 0],
            x_opt['tibialis'][2, 0], x_opt['soleus'][2, 0], x_opt['gastrocnemius'][2, 0],
        ]).flatten()

        all_state = np.concatenate([neuromusculoskeletal_state_num, rooted_variables])

        tendon_force = casadi_function['get_tendon_force'](all_state, mtu_params)
        fiber_force = casadi_function['get_fiber_force'](all_state, mtu_params)
        fiber_passive_force = casadi_function['get_fiber_passive_force'](all_state, mtu_params)
        fiber_active_force = casadi_function['get_fiber_active_force'](all_state, mtu_params)
        ankle_torque = casadi_function['get_joint_moment'](all_state, mtu_params)

        # === diagrammes d'architecture ===
        for axm, name, color in zip(ax_arch, muscle_names, arch_colors):
            xo = x_opt[name]
            draw_muscle_architecture(axm, xo[0, 0], xo[1, 0], xo[2, 0], name, color)

        # === rapport texte (remplace les print) ===
        def g(x, i):  # extraction scalaire sûre
            return float(np.array(x).flatten()[i])

        lines = []
        lines.append('=' * 52)
        lines.append('  PARAMÈTRES MTU COURANTS  [l0m, phi0, f0m, km, lst, kt]')
        lines.append('=' * 52)
        for m, name in enumerate(muscle_names):
            vals = []
            for p in range(n_par):
                v = mtu_params[pidx(p, m)]
                if param_is_angle[p]:
                    vals.append(f'{param_names[p]}={np.rad2deg(v):.2f}°')
                else:
                    vals.append(f'{param_names[p]}={v:.4f}')
            lines.append(f'{name:>13}: ' + '  '.join(vals))

        lines.append('')
        lines.append('=' * 52)
        lines.append('  CONVERGENCE DU SOLVEUR')
        lines.append('=' * 52)
        for name in muscle_names:
            res = results[name]
            lines.append(f'{name:>13}:  convergé={fmt_success(res):<6}  résidu={fmt_residual(res)}')

        lines.append('')
        lines.append('=' * 52)
        lines.append('  DÉTAIL PAR MUSCLE')
        lines.append('=' * 52)
        for m, name in enumerate(muscle_names):
            xo = x_opt[name]
            lines.append(f'--- {name} ---')
            lines.append(f'  MTU length        : {g(mtu_length, m):.4f} m')
            lines.append(f'  moment arm        : {float(mtu_moment_arm[m, -1]):.4f} m')
            lines.append(f'  fiber length      : {g(xo[0, 0], 0):.4f} m')
            lines.append(f'  pennation angle   : {np.rad2deg(g(xo[1, 0], 0)):.4f} deg')
            lines.append(f'  tendon length     : {g(xo[2, 0], 0):.4f} m')
            lines.append(f'  tendon force      : {g(tendon_force, m):.4f} N')
            lines.append(f'  fiber force       : {g(fiber_force, m):.4f} N')
            lines.append(f'    - active        : {g(fiber_active_force, m):.4f} N')
            lines.append(f'    - passive       : {g(fiber_passive_force, m):.4f} N')

        lines.append('')
        lines.append('=' * 52)
        lines.append('  COUPLE ARTICULAIRE')
        lines.append('=' * 52)
        lines.append(f'  ankle moment      : {g(ankle_torque, 0):.4f} N.m')

        set_report('\n'.join(lines))
        fig.canvas.draw_idle()

    # =================================================================
    # CONTRÔLES (fenêtre dédiée « Contrôles ») — trois zones empilées :
    #   A) Cinématique       q1..q6
    #   B) Activation        a1..a3
    #   C) Muscle parameter  l0m..kt (radio muscle + 6 sliders)
    # + bouton Reset
    # IMPORTANT : on attache explicitement les axes à fig_ctrl
    # (sinon ils iraient sur la dernière figure active).
    # =================================================================
    axcolor = 'lightgoldenrodyellow'

    # en-têtes de zone
    fig_ctrl.text(0.06, 0.965, 'A — Cinématique (q1…q6)', fontsize=11, fontweight='bold')
    fig_ctrl.text(0.06, 0.605, 'B — Activation (a1…a3)', fontsize=11, fontweight='bold')
    fig_ctrl.text(0.06, 0.435, 'C — Muscle parameter', fontsize=11, fontweight='bold')

    # ---------------------------------------------------------------
    # ZONE A : cinématique q1..q6
    # ---------------------------------------------------------------
    slider_limits = [(-1, 1), (-1, 1), (-1, 1), (-90, 90), (0, 90), (-40, 40)]
    slider_q = []
    for i, (min_val, max_val) in enumerate(slider_limits):
        # empilés du haut (q1) vers le bas (q6)
        y = 0.915 - i * 0.050
        ax_slider = fig_ctrl.add_axes([0.18, y, 0.70, 0.025], facecolor=axcolor)
        s = Slider(ax_slider, f'q{i + 1}', min_val, max_val, valinit=0.0)
        s.on_changed(update_plot)
        slider_q.append(s)

    # ---------------------------------------------------------------
    # ZONE B : activations a1..a3
    # ---------------------------------------------------------------
    muscle_sliders = []
    for i in range(3):
        y = 0.555 - i * 0.045
        ax_muscle = fig_ctrl.add_axes([0.18, y, 0.70, 0.025], facecolor='mistyrose')
        ms = Slider(ax_muscle, short_names[i], 0.0, 1.0, valinit=0.0, color='red')
        ms.on_changed(update_plot)
        muscle_sliders.append(ms)

    # ---------------------------------------------------------------
    # ZONE C : muscle parameter (radio muscle + 6 sliders)
    # ---------------------------------------------------------------
    def bounds(p, m):
        v0 = muscle_tendon_parameters_num[pidx(p, m)]
        if param_is_angle[p]:                 # phi0 : ±50 %
            lo, hi = v0 * 0.5, v0 * 1.5
        else:                                 # autres : ×0.1 … ×10
            lo, hi = v0 * 0.1, v0 * 10.0
        if lo > hi:
            lo, hi = hi, lo
        if v0 == 0:
            lo, hi = -1.0, 1.0
        return lo, hi

    state = {'muscle': 0}

    # radio de sélection du muscle (gauche)
    ax_radio = fig_ctrl.add_axes([0.06, 0.16, 0.26, 0.24], facecolor=axcolor)
    radio = RadioButtons(ax_radio, muscle_names, active=0)
    ax_radio.set_title('Muscle édité', fontsize=9)

    # 6 sliders MTU (droite)
    mtu_param_sliders = []
    for p in range(n_par):
        y = 0.385 - p * 0.042
        ax_p = fig_ctrl.add_axes([0.50, y, 0.42, 0.022], facecolor='lavender')
        lo, hi = bounds(p, state['muscle'])
        v0 = mtu_params[pidx(p, state['muscle'])]
        label = param_names[p] + ('(°)' if param_is_angle[p] else '')
        disp0 = np.rad2deg(v0) if param_is_angle[p] else v0
        disp_lo = np.rad2deg(lo) if param_is_angle[p] else lo
        disp_hi = np.rad2deg(hi) if param_is_angle[p] else hi
        sp = Slider(ax_p, label, disp_lo, disp_hi, valinit=disp0, color='slateblue')
        mtu_param_sliders.append(sp)

    def on_mtu_change(val=None):
        m = state['muscle']
        for p in range(n_par):
            disp = mtu_param_sliders[p].val
            v = np.deg2rad(disp) if param_is_angle[p] else disp
            mtu_params[pidx(p, m)] = v
        update_plot()

    for sp in mtu_param_sliders:
        sp.on_changed(on_mtu_change)

    def refresh_mtu_sliders(m):
        """recâble bornes + valeurs des 6 sliders sur le muscle m (sans recalcul)."""
        for p in range(n_par):
            lo, hi = bounds(p, m)
            v = mtu_params[pidx(p, m)]
            disp = np.rad2deg(v) if param_is_angle[p] else v
            disp_lo = np.rad2deg(lo) if param_is_angle[p] else lo
            disp_hi = np.rad2deg(hi) if param_is_angle[p] else hi
            sp = mtu_param_sliders[p]
            sp.valmin, sp.valmax = disp_lo, disp_hi
            sp.ax.set_xlim(disp_lo, disp_hi)
            sp.eventson = False
            sp.set_val(disp)
            sp.eventson = True

    def on_muscle_select(label):
        m = muscle_names.index(label)
        state['muscle'] = m
        refresh_mtu_sliders(m)
        fig_ctrl.canvas.draw_idle()

    radio.on_clicked(on_muscle_select)

    # ---------------------------------------------------------------
    # BOUTON RESET
    # ---------------------------------------------------------------
    def reset_all(event=None):
        mtu_params[:] = muscle_tendon_parameters_num
        for s in slider_q + muscle_sliders:
            s.eventson = False
            s.reset()
            s.eventson = True
        for sp in mtu_param_sliders:
            sp.eventson = False
        refresh_mtu_sliders(state['muscle'])
        for sp in mtu_param_sliders:
            sp.eventson = True
        update_plot()

    ax_reset = fig_ctrl.add_axes([0.06, 0.04, 0.26, 0.06])
    btn_reset = Button(ax_reset, 'Reset', color='lightcoral', hovercolor='salmon')
    btn_reset.on_clicked(reset_all)

    update_plot()
    ax.view_init(elev=90, azim=-90)
    plt.show()

def plot_force_length_activation_3d(
    a,
    fiber_length,
    tendon_length,
    muscle_tendon_parameters,
    get_fiber_force_from_fiber_length,
    get_tendon_force_from_tendon_length,
    muscle_idx: int = 1,
    muscle_names: list = None,
    n_grid: int = 50,
    normalize: bool = True,
    elev: int = 15,
    azim: int = -180,
    cmap: str = "cividis",
):
    """
    Représentation 3D force-longueur-activation pour un muscle d'un système
    musculo-tendineux multi-muscles (vecteurs CasADi de taille 3).

    Convention des paramètres (12) :
        [l0m_1, l0m_2, l0m_3,
         phi0_1, phi0_2, phi0_3,
         f0m_1, f0m_2, f0m_3,
         lst_1, lst_2, lst_3]

    Convention des indices muscles : 1, 2, 3 (1-indexé).
    """
    # ---------- coercition en arrays numpy ----------
    a = np.asarray(a, dtype=float).flatten()
    fiber_length = np.asarray(fiber_length, dtype=float).flatten()
    tendon_length = np.asarray(tendon_length, dtype=float).flatten()
    params = np.asarray(muscle_tendon_parameters, dtype=float).flatten()

    n_muscles = 3
    assert a.size == n_muscles
    assert fiber_length.size == n_muscles
    assert tendon_length.size == n_muscles
    assert params.size == 12
    assert muscle_idx in (1, 2, 3), "muscle_idx doit être 1, 2 ou 3"

    if muscle_names is None:
        muscle_names = [f"Muscle {i}" for i in (1, 2, 3)]

    # ---------- extraction des paramètres du muscle d'intérêt ----------
    # Paramètres rangés par type : [l0m×3, phi0×3, f0m×3, lst×3]
    p_offset = muscle_idx - 1
    l0m = float(params[p_offset + 0])
    f0m = float(params[p_offset + 6])
    lst = float(params[p_offset + 9])

    # index 0-based pour accéder aux vecteurs d'état (a, fiber_length, tendon_length)
    m = muscle_idx - 1

    # ---------- helpers d'évaluation des fonctions CasADi vectorielles ----------
    def eval_fiber_force(a_scalar, l_scalar):
        a_vec = a.copy()
        l_vec = fiber_length.copy()
        a_vec[m] = a_scalar
        l_vec[m] = l_scalar
        out = np.array(
            get_fiber_force_from_fiber_length(a_vec, l_vec, params)
        ).flatten()
        return float(out[m])

    def eval_tendon_force(lt_scalar):
        lt_vec = tendon_length.copy()
        lt_vec[m] = lt_scalar
        out = np.array(
            get_tendon_force_from_tendon_length(lt_vec, params)
        ).flatten()
        return float(out[m])

    # ---------- grille muscle ----------
    fl_range = np.linspace(0.4 * l0m, 1.8 * l0m, n_grid)
    a_range = np.linspace(0.0, 1.0, n_grid)
    A_grid, L_grid = np.meshgrid(a_range, fl_range)

    passive_force = np.zeros_like(L_grid)
    active_force = np.zeros_like(L_grid)
    total_force = np.zeros_like(L_grid)

    for i in range(L_grid.shape[0]):
        for j in range(L_grid.shape[1]):
            l = L_grid[i, j]
            act = A_grid[i, j]
            try:
                f_pass = eval_fiber_force(0.0, l)
                f_tot = eval_fiber_force(act, l)
                passive_force[i, j] = f_pass
                total_force[i, j] = f_tot
                active_force[i, j] = f_tot - f_pass
            except Exception as e:
                print(f"[fiber] erreur i={i}, j={j}, l={l}, a={act}: {e}")
                passive_force[i, j] = np.nan
                active_force[i, j] = np.nan
                total_force[i, j] = np.nan

    # ---------- point de fonctionnement courant ----------
    a_pt = a[m]
    fl_pt_raw = fiber_length[m]
    f_pass_pt = eval_fiber_force(0.0, fl_pt_raw)
    f_tot_pt = eval_fiber_force(a_pt, fl_pt_raw)
    f_act_pt = f_tot_pt - f_pass_pt

    # ---------- courbe tendon ----------
    lt_range = np.linspace(0.97 * lst, 1.10 * lst, n_grid * 4)
    tendon_curve = np.array([eval_tendon_force(lt) for lt in lt_range])
    lt_pt_raw = tendon_length[m]
    f_tendon_pt = eval_tendon_force(lt_pt_raw)

    # ---------- normalisation ----------
    if normalize:
        L_plot = L_grid / l0m
        fl_pt = fl_pt_raw / l0m
        passive_plot = passive_force / f0m
        active_plot = active_force / f0m
        total_plot = total_force / f0m
        f_pass_pt_plot = f_pass_pt / f0m
        f_act_pt_plot = f_act_pt / f0m
        f_tot_pt_plot = f_tot_pt / f0m

        lt_plot = lt_range / lst
        lt_pt_plot = lt_pt_raw / lst
        tendon_curve_plot = tendon_curve / f0m
        f_tendon_pt_plot = f_tendon_pt / f0m

        ylabel_l = r"$\tilde{l}_f = l_f / l_0^m$"
        zlabel_F = r"$\tilde{F} = F / F_0^m$"
        xlabel_lt = r"$\tilde{l}_t = l_t / l_{st}$"
    else:
        L_plot = L_grid
        fl_pt = fl_pt_raw
        passive_plot, active_plot, total_plot = passive_force, active_force, total_force
        f_pass_pt_plot, f_act_pt_plot, f_tot_pt_plot = f_pass_pt, f_act_pt, f_tot_pt

        lt_plot = lt_range
        lt_pt_plot = lt_pt_raw
        tendon_curve_plot = tendon_curve
        f_tendon_pt_plot = f_tendon_pt

        ylabel_l = "Longueur de fibre (m)"
        zlabel_F = "Force (N)"
        xlabel_lt = "Longueur de tendon (m)"

    # ---------- figure ----------
    fig = plt.figure(figsize=(20, 5.2))
    fig.suptitle(
        f"Force – longueur – activation : {muscle_names[m]}",
        fontsize=13, fontweight="bold", y=1.02,
    )

    titles = ["Force passive", "Force active", "Force totale (passive + active)"]
    surfaces = [passive_plot, active_plot, total_plot]
    points_z = [f_pass_pt_plot, f_act_pt_plot, f_tot_pt_plot]
    axes_3d = []

    for k, (title, F_grid, z_pt) in enumerate(zip(titles, surfaces, points_z)):
        ax = fig.add_subplot(1, 4, k + 1, projection="3d")
        surf = ax.plot_surface(
            A_grid, L_plot, F_grid,
            cmap=cmap, edgecolor="none", alpha=0.92, antialiased=True,
        )
        ax.scatter(
            [a_pt], [fl_pt], [z_pt],
            color="k", s=80, marker="*", depthshade=False,
            edgecolor="white", linewidth=0.8, zorder=10,
        )
        zmin = np.nanmin(F_grid)
        ax.plot([a_pt, a_pt], [fl_pt, fl_pt], [zmin, z_pt],
                color="k", lw=0.6, ls=":")

        ax.set_title(title, fontsize=11, pad=8)
        ax.set_xlabel("Activation $a$", fontsize=9, labelpad=4)
        ax.set_ylabel(ylabel_l, fontsize=9, labelpad=4)
        ax.set_zlabel(zlabel_F, fontsize=9, labelpad=4)
        ax.tick_params(labelsize=8)
        ax.view_init(elev=elev, azim=azim)
        ax.invert_yaxis()
        fig.colorbar(surf, ax=ax, shrink=0.55, aspect=14, pad=0.08)
        axes_3d.append(ax)

    # ---------- panneau tendon ----------
    ax_t = fig.add_subplot(1, 4, 4)
    ax_t.fill_between(lt_plot, 0, tendon_curve_plot, color="#BA7517", alpha=0.15)
    ax_t.plot(lt_plot, tendon_curve_plot, color="#BA7517", lw=2.4, label="Force tendon")
    ax_t.scatter([lt_pt_plot], [f_tendon_pt_plot], color="k", marker="*",
                 s=110, zorder=5, label="Point courant")
    ax_t.set_title("Force tendon", fontsize=11, pad=8)
    ax_t.set_xlabel(xlabel_lt, fontsize=10)
    ax_t.set_ylabel(zlabel_F, fontsize=10)
    ax_t.tick_params(labelsize=9)
    ax_t.grid(True, alpha=0.3, lw=0.6)
    ax_t.spines[["top", "right"]].set_visible(False)
    ax_t.legend(fontsize=8.5, frameon=False, loc="upper left")

    plt.tight_layout()
    return fig, (*axes_3d, ax_t)

def plot_all_muscles_3d(
    a,
    fiber_length,
    tendon_length,
    muscle_tendon_parameters,
    get_fiber_force_from_fiber_length,
    get_tendon_force_from_tendon_length,
    muscle_names: list = None,
    **kwargs,
):
    """Boucle sur les 3 muscles (indices 1, 2, 3)."""
    if muscle_names is None:
        muscle_names = [f"Muscle {i}" for i in (1, 2, 3)]

    figures = []
    for idx in (1, 2, 3):
        fig, _ = plot_force_length_activation_3d(
            a=a,
            fiber_length=fiber_length,
            tendon_length=tendon_length,
            muscle_tendon_parameters=muscle_tendon_parameters,
            get_fiber_force_from_fiber_length=get_fiber_force_from_fiber_length,
            get_tendon_force_from_tendon_length=get_tendon_force_from_tendon_length,
            muscle_idx=idx,
            muscle_names=muscle_names,
            **kwargs,
        )
        figures.append(fig)
    return figures



def force_length_sliders(
    param_config: dict,
    muscle_tendon_parameters_num,
    a,
    fiber_length,
    tendon_length,
    get_fiber_force_from_fiber_length,
    get_tendon_force_from_tendon_length,
    muscle_idx: int = 1,
    muscle_names: list = None,
    mtu_params=("l0m", "phi0", "f0m", "km", "lst", "kt"),
    n_muscles: int = 3,
    param_layout: str = "by_type",   # "by_type" (idx = t*n_muscles + m) ou "by_muscle"
    n_grid: int = 200,
    normalize: bool = True,
    slider_ranges: dict = None,
):
    """
    Paramètres
    ----------
    param_config : {nom: 'sym' | 'fixed'}
        'sym'   -> slider (valeur modifiable, passée aux équations).
        'fixed' -> pas de slider (valeur figée = celle du vecteur 18).
    muscle_tendon_parameters_num : array (18,)
        Vecteur complet des paramètres (6 types × 3 muscles).
    a, fiber_length, tendon_length : array (3,)
        États initiaux des 3 muscles.
    muscle_idx : 1, 2 ou 3
        Muscle visualisé.
    param_layout : 'by_type' | 'by_muscle'
        Rangement du vecteur 18 (voir pidx ci-dessous).
    slider_ranges : {nom: (min, max, step)}  surcharge des bornes par défaut.

    Retour
    ------
    dict des widgets (à GARDER dans une variable, sinon ils sont désactivés
    par le ramasse-miettes).
    """
    mtu_params = list(mtu_params)
    for need in ("l0m", "f0m", "lst"):
        assert need in mtu_params, f"mtu_params doit contenir '{need}'"
    assert muscle_idx in range(1, n_muscles + 1)

    params_base = np.asarray(muscle_tendon_parameters_num, float).flatten().copy()
    a_base = np.asarray(a, float).flatten().copy()
    fl_base = np.asarray(fiber_length, float).flatten().copy()
    tl_base = np.asarray(tendon_length, float).flatten().copy()
    assert params_base.size == n_muscles * len(mtu_params), (
        f"params attendu de taille {n_muscles*len(mtu_params)}, reçu {params_base.size}"
    )
    for v, nm in ((a_base, "a"), (fl_base, "fiber_length"), (tl_base, "tendon_length")):
        assert v.size == n_muscles, f"{nm} doit être de taille {n_muscles}"

    if muscle_names is None:
        muscle_names = [f"Muscle {i}" for i in range(1, n_muscles + 1)]

    state = {"m": muscle_idx - 1, "normalize": bool(normalize)}
    sym_names = [p for p in mtu_params if param_config.get(p) == "sym"]

    # ---------- index d'un paramètre (type `name`, muscle 0-based `m0`) ----------
    def pidx(name, m0):
        t = mtu_params.index(name)
        if param_layout == "by_type":
            return t * n_muscles + m0
        elif param_layout == "by_muscle":
            return m0 * len(mtu_params) + t
        raise ValueError("param_layout doit être 'by_type' ou 'by_muscle'")

    # ---------- valeur courante d'un paramètre du muscle sélectionné ----------
    def cur(name):
        if param_config.get(name) == "sym":
            return float(sym_sliders[name].val)
        return float(params_base[pidx(name, state["m"])])

    # ---------- vecteur params (18) reconstruit pour évaluation ----------
    def build_params():
        p = params_base.copy()
        for name in sym_names:
            p[pidx(name, state["m"])] = sym_sliders[name].val
        return p

    # ---------- évaluateurs scalaires (sur le muscle sélectionné) ----------
    def eval_fiber(a_val, l_val, p):
        m = state["m"]
        av = a_base.copy(); av[m] = a_val
        lv = fl_base.copy(); lv[m] = l_val
        out = np.array(get_fiber_force_from_fiber_length(av, lv, p)).flatten()
        return float(out[m])

    def eval_tendon(lt_val, p):
        m = state["m"]
        lv = tl_base.copy(); lv[m] = lt_val
        out = np.array(get_tendon_force_from_tendon_length(lv, p)).flatten()
        return float(out[m])

    # ---------- bornes par défaut des sliders sym (valides pour les 3 muscles) ----------
    def default_range(name):
        if slider_ranges and name in slider_ranges:
            return slider_ranges[name]
        if name == "phi0":
            return 0.0, np.pi / 2, 0.005
        vals = np.array([params_base[pidx(name, mm)] for mm in range(n_muscles)])
        vmin, vmax = float(vals.min()), float(vals.max())
        if vmin > 0:
            lo, hi = 0.1 * vmin, 3 * vmax
        elif vmax < 0:
            lo, hi = 3 * vmin, 0.1 * vmax
        else:
            span = max(abs(vmin), abs(vmax), 1.0)
            lo, hi = vmin - span, vmax + span
        return lo, hi, (hi - lo) / 200.0

    # ==================== figure & widgets ====================
    fig = plt.figure(figsize=(14, 9))
    ax_m = fig.add_axes([0.06, 0.62, 0.40, 0.32])
    ax_t = fig.add_axes([0.55, 0.62, 0.40, 0.32])
    txt_force = fig.text(0.5, 0.975, "", ha="center", fontsize=10)

    fig.text(0.10, 0.555, "Paramètres MTU (muscle sélectionné)",
             fontweight="bold", fontsize=10)
    fig.text(0.60, 0.555, "Point de fonctionnement", fontweight="bold", fontsize=10)

    def make_slider(rect, label, vmin, vmax, vinit, fmt, step=None):
        ax = fig.add_axes(rect)
        return Slider(ax, label, vmin, vmax, valinit=vinit, valfmt=fmt, valstep=step)

    # --- sliders paramètres sym ---
    m0 = state["m"]
    sym_sliders = {}
    for k, name in enumerate(sym_names):
        lo, hi, step = default_range(name)
        v0 = float(params_base[pidx(name, m0)])
        sym_sliders[name] = make_slider(
            [0.12, 0.50 - 0.05 * k, 0.30, 0.022], name, lo, hi, v0, "%.4g"
        )

    # --- sliders point de fonctionnement (longueurs en NORMALISÉ) ---
    l0m0 = float(params_base[pidx("l0m", m0)])
    lst0 = float(params_base[pidx("lst", m0)])
    s_a = make_slider([0.62, 0.50, 0.26, 0.022], "a (activ.)", 0.0, 1.0,
                      float(a_base[m0]), "%.2f")
    s_lf = make_slider([0.62, 0.45, 0.26, 0.022], r"$\tilde{l}_f$", 0.5, 1.6,
                       float(fl_base[m0]) / l0m0, "%.3f")
    s_lt = make_slider([0.62, 0.40, 0.26, 0.022], r"$\tilde{l}_t$", 0.97, 1.10,
                       float(tl_base[m0]) / lst0, "%.4f")

    # --- RadioButtons muscle ---
    fig.text(0.62, 0.37, "Muscle", fontweight="bold", fontsize=10)
    ax_radio_m = fig.add_axes([0.62, 0.20, 0.15, 0.16])
    radio_m = RadioButtons(ax_radio_m, muscle_names, active=m0)

    # --- RadioButtons affichage ---
    fig.text(0.80, 0.31, "Affichage", fontweight="bold", fontsize=10)
    ax_radio_n = fig.add_axes([0.80, 0.20, 0.15, 0.10])
    radio_n = RadioButtons(ax_radio_n, ("Normalisé", "Physique"),
                           active=0 if normalize else 1)

    # --- Button reset ---
    ax_btn = fig.add_axes([0.80, 0.07, 0.10, 0.05])
    btn_reset = Button(ax_btn, "Reset")

    # ---------- artistes (créés une fois) ----------
    (l_pass,) = ax_m.plot([], [], "--", color="#1f77b4", lw=2, label="Passive")
    (l_act,) = ax_m.plot([], [], "-", color="#d62728", lw=2, label="Active")
    (l_tot,) = ax_m.plot([], [], "-", color="#2ca02c", lw=2.4, label="Totale")
    sc_m = ax_m.scatter([], [], color="k", marker="*", s=90, zorder=5)
    vline_m = ax_m.axvline(np.nan, color="k", ls=":", lw=0.7)
    ax_m.grid(True, alpha=0.3)
    ax_m.spines[["top", "right"]].set_visible(False)
    ax_m.legend(fontsize=9, frameon=False, loc="upper left")

    (l_ten,) = ax_t.plot([], [], color="#BA7517", lw=2.4, label="Tendon")
    sc_t = ax_t.scatter([], [], color="k", marker="*", s=110, zorder=5)
    vline_t = ax_t.axvline(np.nan, color="k", ls=":", lw=0.7)
    ax_t.grid(True, alpha=0.3)
    ax_t.spines[["top", "right"]].set_visible(False)
    ax_t.legend(fontsize=9, frameon=False, loc="upper left")
    holder = {"fill": None}

    # ==================== redraw ====================
    def redraw(_=None):
        m = state["m"]
        norm = state["normalize"]
        p = build_params()
        l0m, f0m, lst = cur("l0m"), cur("f0m"), cur("lst")
        a_val, lf_n, lt_n = s_a.val, s_lf.val, s_lt.val
        lf_phys, lt_phys = lf_n * l0m, lt_n * lst

        lf_grid_n = np.linspace(0.5, 1.6, n_grid)
        lf_grid = lf_grid_n * l0m
        passive = np.array([eval_fiber(0.0, l, p) for l in lf_grid])
        total = np.array([eval_fiber(a_val, l, p) for l in lf_grid])
        active = total - passive

        lt_grid_n = np.linspace(0.97, 1.10, n_grid)
        lt_grid = lt_grid_n * lst
        tendon = np.array([eval_tendon(l, p) for l in lt_grid])

        f_pass = eval_fiber(0.0, lf_phys, p)
        f_tot = eval_fiber(a_val, lf_phys, p)
        f_act = f_tot - f_pass
        f_ten = eval_tendon(lt_phys, p)

        if norm:
            xf, xt = lf_grid_n, lt_grid_n
            yp, ya, yt = passive / f0m, active / f0m, total / f0m
            yten = tendon / f0m
            xf_op, xt_op = lf_n, lt_n
            yp_op, ya_op, yt_op, yten_op = f_pass / f0m, f_act / f0m, f_tot / f0m, f_ten / f0m
            ax_m.set_xlabel(r"$\tilde{l}_f = l_f/l_0^m$", fontsize=10)
            ax_t.set_xlabel(r"$\tilde{l}_t = l_t/l_{st}$", fontsize=10)
            ax_m.set_ylabel(r"$\tilde{F} = F/F_0^m$", fontsize=10)
            ax_t.set_ylabel(r"$\tilde{F} = F/F_0^m$", fontsize=10)
        else:
            xf, xt = lf_grid, lt_grid
            yp, ya, yt = passive, active, total
            yten = tendon
            xf_op, xt_op = lf_phys, lt_phys
            yp_op, ya_op, yt_op, yten_op = f_pass, f_act, f_tot, f_ten
            ax_m.set_xlabel("Longueur de fibre (m)", fontsize=10)
            ax_t.set_xlabel("Longueur de tendon (m)", fontsize=10)
            ax_m.set_ylabel("Force (N)", fontsize=10)
            ax_t.set_ylabel("Force (N)", fontsize=10)

        l_pass.set_data(xf, yp)
        l_act.set_data(xf, ya)
        l_tot.set_data(xf, yt)
        sc_m.set_offsets(np.c_[[xf_op] * 3, [yp_op, ya_op, yt_op]])
        vline_m.set_xdata([xf_op, xf_op])

        l_ten.set_data(xt, yten)
        sc_t.set_offsets(np.c_[[xt_op], [yten_op]])
        vline_t.set_xdata([xt_op, xt_op])
        if holder["fill"] is not None:
            holder["fill"].remove()
        holder["fill"] = ax_t.fill_between(xt, 0, yten, color="#BA7517", alpha=0.15)

        # --- bornes (calculées à la main : scatter/fill non pris par autoscale) ---
        def set_lims(ax, xv, ys):
            xmin, xmax = float(np.min(xv)), float(np.max(xv))
            ymin = float(np.min([np.min(y) for y in ys]))
            ymax = float(np.max([np.max(y) for y in ys]))
            ylo = min(0.0, ymin)
            pad = 0.05 * (ymax - ylo) if ymax > ylo else 1.0
            ax.set_xlim(xmin, xmax)
            ax.set_ylim(ylo - pad, ymax + pad)

        set_lims(ax_m, xf, [yp, ya, yt])
        set_lims(ax_t, xt, [yten, [0.0]])

        ax_m.set_title(f"{muscle_names[m]} — fibre (a = {a_val:.2f})", fontsize=11)
        ax_t.set_title(f"{muscle_names[m]} — tendon", fontsize=11)
        txt_force.set_text(
            f"F_passive = {f_pass:.1f} N    F_active = {f_act:.1f} N    "
            f"F_totale = {f_tot:.1f} N    F_tendon = {f_ten:.1f} N"
        )
        fig.canvas.draw_idle()

    # ==================== callbacks ====================
    for s in (*sym_sliders.values(), s_a, s_lf, s_lt):
        s.on_changed(redraw)

    def on_muscle(label):
        state["m"] = muscle_names.index(label)
        m = state["m"]
        l0m_m = float(params_base[pidx("l0m", m)])
        lst_m = float(params_base[pidx("lst", m)])
        for name in sym_names:                       # déclenche redraw via on_changed
            sym_sliders[name].set_val(float(params_base[pidx(name, m)]))
        s_a.set_val(float(a_base[m]))
        s_lf.set_val(float(fl_base[m]) / l0m_m)
        s_lt.set_val(float(tl_base[m]) / lst_m)
        redraw()

    def on_norm(label):
        state["normalize"] = (label == "Normalisé")
        redraw()

    def on_reset(_event):
        m = state["m"]
        l0m_m = float(params_base[pidx("l0m", m)])
        lst_m = float(params_base[pidx("lst", m)])
        for name in sym_names:
            sym_sliders[name].set_val(float(params_base[pidx(name, m)]))
        s_a.set_val(float(a_base[m]))
        s_lf.set_val(float(fl_base[m]) / l0m_m)
        s_lt.set_val(float(tl_base[m]) / lst_m)

    radio_m.on_clicked(on_muscle)
    radio_n.on_clicked(on_norm)
    btn_reset.on_clicked(on_reset)

    redraw()
    plt.show()

    return {
        "fig": fig,
        "sym_sliders": sym_sliders,
        "s_a": s_a, "s_lf": s_lf, "s_lt": s_lt,
        "radio_muscle": radio_m, "radio_norm": radio_n, "btn_reset": btn_reset,
    }





def interactive_calibration_figure(
        data,
        skeleton_num,
        casadi_function,
        param_index,
        param_initial,
        param_lower=None,
        param_upper=None,
        muscles=('ta', 'sol', 'gast'),
        solve_all_muscles=None,
):
    """
    Figure interactive : sliders sur les paramètres muscle-tendon, et tracés
    en temps réel de mesuré / simulé / différence + coûts (torque, fibre, total)
    dans les deux formulations (pondérée et non pondérée).

    IMPORTANT — cohérence avec ton modèle
    -------------------------------------
    Le vecteur `muscle_tendon_parameters` attendu par les fonctions CasADi
    dépend du `param_config` utilisé dans `get_model_equation`. Son layout est
    PARAM-MAJOR : un bloc [ta, sol, gast] par paramètre 'sym', dans l'ordre de
    `param_index`. Pour piloter `kt` avec un slider, le `casadi_function` fourni
    DOIT avoir été construit avec kt en 'sym', sinon kt est codé en dur dans les
    équations et le slider n'aura aucun effet.

    Reconstruction recommandée AVANT d'appeler cette fonction :

        param_config = {'l0m':'sym','phi0':'sym','f0m':'sym',
                        'lst':'sym','kt':'sym','km':'fixed'}
        casadi_function, mtu_sym, _ = get_model_equation(
            param_config=param_config, fixed_values={'km':4.0})
        initial_guess, ub, lb, param_index = get_initial_guess(
            mtu_ref, data, param_config, strategy='measured')

    Parameters
    ----------
    data : np.ndarray, shape (15, n_trials)
        Données mesurées de référence (mêmes unités SI que dans ton pipeline :
        torque N.m, angles rad, longueurs m). Lignes :
          [0]=torque, [1:3]=q(knee,ankle), [3:6]=activations,
          [6:9]=fiber length, [9:12]=pennation, [12:15]=tendon length.
    skeleton_num : array-like
        Géométrie squelette.
    casadi_function : dict
        Fonctions CasADi (cf. get_model_equation). Doit contenir
        'get_mtu_length', 'get_joint_moment', et le rootfinder via
        solve_all_muscles.
    param_index : dict {param_name: {muscle: idx}}
        Mapping paramètre/muscle -> index dans le vecteur optimisé. Donne
        l'ordre exact du vecteur muscle_tendon_parameters.
    param_initial : np.ndarray, shape (n_param,)
        Valeurs initiales des paramètres (= position de départ des sliders),
        dans l'ordre de param_index.
    param_lower, param_upper : np.ndarray, shape (n_param,), optional
        Bornes des sliders. Défaut : [0.5*init, 1.5*init] (et garde-fous
        physiologiques pour phi0).
    muscles : tuple of str
        Noms des muscles, ordre [ta, sol, gast].
    solve_all_muscles : callable
        Ta fonction solve_all_muscles(a, mtu_length, mtu_params, casadi_function).
        Obligatoire (passée pour éviter les imports circulaires).

    Returns
    -------
    fig : matplotlib Figure (gardée vivante par plt.show()).
    """
    if solve_all_muscles is None:
        raise ValueError("Passe ta fonction solve_all_muscles en argument.")

    data = np.asarray(data, dtype=float)
    assert data.shape[0] == 15, f"data doit être (15, n_trials), reçu {data.shape}"
    n_trials = data.shape[1]

    param_initial = np.asarray(param_initial, dtype=float)
    n_param = param_initial.shape[0]

    # --- ordre à plat des (param, muscle) d'après param_index ---
    flat_keys = [None] * n_param
    for p, mp in param_index.items():
        for m, idx in mp.items():
            flat_keys[idx] = (p, m)
    assert all(k is not None for k in flat_keys), \
        "param_index incohérent avec param_initial (indices manquants)."

    # --- bornes des sliders ---
    if param_lower is None or param_upper is None:
        lo = np.empty(n_param)
        hi = np.empty(n_param)
        for idx, (p, m) in enumerate(flat_keys):
            v = param_initial[idx]
            if p == 'phi0':                       # pennation : ±15° autour, dans [1°, 45°]
                lo[idx] = max(v - np.deg2rad(15), np.deg2rad(1))
                hi[idx] = min(v + np.deg2rad(15), np.deg2rad(45))
            else:
                lo[idx] = 0.5 * v
                hi[idx] = 1.5 * v
        param_lower = lo if param_lower is None else np.asarray(param_lower, float)
        param_upper = hi if param_upper is None else np.asarray(param_upper, float)

    # --- poids du coût pondéré (mêmes valeurs que optimization_nlp) ---
    sigma_torque = 10         # N.m
    sigma_fiber = 0.01         # m

    # --- mesures de référence (lignes utiles) ---
    meas_torque = data[0, :]
    meas_fiber = data[6:9, :]   # (3, n_trials)  ta, sol, gast

    trials = np.arange(n_trials)

    # ============ Cœur : simulation forward pour un vecteur params ============ #
    def simulate(mtu_params):
        """Renvoie torque simulé (n_trials,) et fiber length simulée (3,n_trials)."""
        sim_torque = np.full(n_trials, np.nan)
        sim_fiber = np.full((3, n_trials), np.nan)

        for t in range(n_trials):
            col = data[:, t]
            a_trial = col[3:6]
            q_num = [0, 0, 0, 0, col[1], col[2]]
            musc_states = q_num + list(skeleton_num)
            neuromusc = np.concatenate([a_trial, musc_states])

            mtu_length = casadi_function['get_mtu_length'](musc_states)

            try:
                results = solve_all_muscles(a_trial, mtu_length, mtu_params,
                                            casadi_function)
            except Exception:
                continue

            x_ta = results['ta']['x_opt'] if 'ta' in results else results['tibialis']['x_opt']
            x_sol = results['sol']['x_opt'] if 'sol' in results else results['soleus']['x_opt']
            x_gast = results['gast']['x_opt'] if 'gast' in results else results['gastrocnemius']['x_opt']

            fl = np.array([x_ta[0, 0], x_sol[0, 0], x_gast[0, 0]]).flatten()
            penn = np.array([x_ta[1, 0], x_sol[1, 0], x_gast[1, 0]]).flatten()
            tl = np.array([x_ta[2, 0], x_sol[2, 0], x_gast[2, 0]]).flatten()

            rooted = np.concatenate([fl, penn, tl])
            all_state = np.concatenate([neuromusc, rooted])

            sim_torque[t] = float(casadi_function['get_joint_moment'](all_state, mtu_params))
            sim_fiber[:, t] = fl

        return sim_torque, sim_fiber

    # ============ Coûts ============ #
    def compute_costs(sim_torque, sim_fiber):
        e_torque = meas_torque - sim_torque              # (n_trials,)
        e_fiber = meas_fiber - sim_fiber                 # (3, n_trials)

        # masque NaN (trials non résolus)
        ok_t = ~np.isnan(e_torque)
        ok_f = ~np.isnan(e_fiber)
        n_ok = max(int(ok_t.sum()), 1)


        # --- pondéré (optimization_nlp) : normalisé par n_element*n_trials ---
        n_element = 4  # 1 torque + 3 fibres
        n_comp = n_element * n_trials
        Jw_torque = np.nansum(e_torque[ok_t] ** 2) / sigma_torque ** 2 / n_comp
        Jw_fiber = np.nansum(e_fiber[ok_f] ** 2) / sigma_fiber ** 2 / n_comp
        Jw_total = Jw_torque + Jw_fiber



        # --- non pondéré (somme des carrés, unités SI) ---
        Ju_torque = np.nansum(e_torque[ok_t] ** 2)
        Ju_fiber = np.nansum(e_fiber[ok_f] ** 2)
        Ju_total = Ju_torque + Ju_fiber

        return {
            'w': (Jw_torque, Jw_fiber, Jw_total),
            'u': (Ju_torque, Ju_fiber, Ju_total),
            'e_torque': e_torque, 'e_fiber': e_fiber, 'n_ok': n_ok,
        }

    # ============ Mise en page de la figure ============ #
    fig = plt.figure(figsize=(15, 9))
    fig.suptitle("Calibration interactive muscle-tendon — mesuré vs simulé",
                 fontsize=13, fontweight='bold')

    # zone de tracé à gauche, sliders à droite
    gs = fig.add_gridspec(3, 2, left=0.06, right=0.62, top=0.92, bottom=0.08,
                          hspace=0.45, wspace=0.28)

    ax_torque = fig.add_subplot(gs[0, 0])
    ax_torque_diff = fig.add_subplot(gs[0, 1])
    ax_fiber = fig.add_subplot(gs[1, 0])
    ax_fiber_diff = fig.add_subplot(gs[1, 1])
    ax_cost_w = fig.add_subplot(gs[2, 0])
    ax_cost_u = fig.add_subplot(gs[2, 1])

    colors = {'ta': 'tab:blue', 'sol': 'tab:green', 'gast': 'tab:red'}

    # --- état initial ---
    sim_torque, sim_fiber = simulate(param_initial)
    C = compute_costs(sim_torque, sim_fiber)

    # ---- (1) torque mesuré vs simulé ----
    ax_torque.set_title("Torque cheville")
    l_meas_t, = ax_torque.plot(trials, meas_torque, 'k.-', lw=1, ms=4, label='mesuré')
    l_sim_t, = ax_torque.plot(trials, sim_torque, 'm.-', lw=1, ms=4, label='simulé')
    ax_torque.set_ylabel("N.m")
    ax_torque.set_xlabel("trial")
    ax_torque.legend(fontsize=8, loc='best')
    ax_torque.grid(alpha=0.3)

    # ---- (2) différence torque ----
    ax_torque_diff.set_title("Différence torque (mesuré - simulé)")
    l_diff_t, = ax_torque_diff.plot(trials, C['e_torque'], 'm.-', lw=1, ms=4)
    ax_torque_diff.axhline(0, color='k', lw=0.6)
    ax_torque_diff.set_ylabel("N.m")
    ax_torque_diff.set_xlabel("trial")
    ax_torque_diff.grid(alpha=0.3)

    # ---- (3) fiber length mesurée vs simulée ----
    ax_fiber.set_title("Fiber length (— mesuré, -- simulé)")
    l_meas_f, l_sim_f = {}, {}
    for i, m in enumerate(muscles):
        l_meas_f[m], = ax_fiber.plot(trials, meas_fiber[i] * 100,
                                     color=colors[m], lw=1.2, label=f'{m} mes.')
        l_sim_f[m], = ax_fiber.plot(trials, sim_fiber[i] * 100,
                                    color=colors[m], lw=1.2, ls='--')
    ax_fiber.set_ylabel("cm")
    ax_fiber.set_xlabel("trial")
    ax_fiber.legend(fontsize=7, loc='best', ncol=3)
    ax_fiber.grid(alpha=0.3)

    # ---- (4) différence fiber length ----
    ax_fiber_diff.set_title("Différence fiber length (mesuré - simulé)")
    l_diff_f = {}
    for i, m in enumerate(muscles):
        l_diff_f[m], = ax_fiber_diff.plot(trials, C['e_fiber'][i] * 100,
                                          color=colors[m], lw=1.2, label=m)
    ax_fiber_diff.axhline(0, color='k', lw=0.6)
    ax_fiber_diff.set_ylabel("cm")
    ax_fiber_diff.set_xlabel("trial")
    ax_fiber_diff.legend(fontsize=7, loc='best', ncol=3)
    ax_fiber_diff.grid(alpha=0.3)

    # ---- (5) coût pondéré ----
    ax_cost_w.set_title("Coût pondéré (optimization_nlp)")
    bars_w = ax_cost_w.bar(['torque', 'fibre', 'total'], C['w'],
                           color=['m', 'tab:cyan', 'tab:gray'])
    ax_cost_w.set_ylabel("J (norm.)")
    ax_cost_w.grid(alpha=0.3, axis='y')
    txt_w = [ax_cost_w.text(b.get_x() + b.get_width() / 2, b.get_height(),
                            f"{v:.3e}", ha='center', va='bottom', fontsize=7)
             for b, v in zip(bars_w, C['w'])]

    # ---- (6) coût non pondéré ----
    ax_cost_u.set_title("Coût non pondéré (Σ carrés, SI)")
    bars_u = ax_cost_u.bar(['torque', 'fibre', 'total'], C['u'],
                           color=['m', 'tab:cyan', 'tab:gray'])
    ax_cost_u.set_ylabel("Σ e²")
    ax_cost_u.grid(alpha=0.3, axis='y')
    txt_u = [ax_cost_u.text(b.get_x() + b.get_width() / 2, b.get_height(),
                            f"{v:.3e}", ha='center', va='bottom', fontsize=7)
             for b, v in zip(bars_u, C['u'])]

    # ============ Sliders (à droite) ============ #
    sliders = []
    slider_left = 0.72
    slider_w = 0.22
    slider_h = 0.018
    y0 = 0.92
    dy = 0.026

    def _fmt_label(p, m):
        if p == 'phi0':
            return f"{p}_{m} [°]"
        if p in ('l0m', 'lst'):
            return f"{p}_{m} [cm]"
        if p == 'f0m':
            return f"{p}_{m} [N]"
        return f"{p}_{m}"

    def _to_display(p, val):
        if p == 'phi0':
            return np.rad2deg(val)
        if p in ('l0m', 'lst'):
            return val * 100
        return val

    def _from_display(p, val):
        if p == 'phi0':
            return np.deg2rad(val)
        if p in ('l0m', 'lst'):
            return val / 100
        return val

    fig.text(slider_left, y0 + 0.03, "Paramètres muscle-tendon",
             fontsize=11, fontweight='bold')

    y = y0
    for idx, (p, m) in enumerate(flat_keys):
        ax_s = fig.add_axes([slider_left, y, slider_w, slider_h])
        s = Slider(
            ax_s,
            _fmt_label(p, m),
            _to_display(p, param_lower[idx]),
            _to_display(p, param_upper[idx]),
            valinit=_to_display(p, param_initial[idx]),
            valfmt="%.3f",
        )
        s.label.set_fontsize(7)
        s.valtext.set_fontsize(7)
        sliders.append(s)
        y -= dy

    # bouton reset
    ax_reset = fig.add_axes([slider_left, y - 0.02, 0.08, 0.03])
    btn_reset = Button(ax_reset, 'Reset')

    # ============ Callback de mise à jour ============ #
    def update(_=None):
        mtu_params = np.array([
            _from_display(flat_keys[i][0], sliders[i].val)
            for i in range(n_param)
        ])

        sim_torque, sim_fiber = simulate(mtu_params)
        C = compute_costs(sim_torque, sim_fiber)

        # torque
        l_sim_t.set_ydata(sim_torque)
        l_diff_t.set_ydata(C['e_torque'])
        ax_torque.relim(); ax_torque.autoscale_view()
        ax_torque_diff.relim(); ax_torque_diff.autoscale_view()

        # fiber
        for i, m in enumerate(muscles):
            l_sim_f[m].set_ydata(sim_fiber[i] * 100)
            l_diff_f[m].set_ydata(C['e_fiber'][i] * 100)
        ax_fiber.relim(); ax_fiber.autoscale_view()
        ax_fiber_diff.relim(); ax_fiber_diff.autoscale_view()

        # coûts pondérés
        for b, v, t in zip(bars_w, C['w'], txt_w):
            b.set_height(v)
            t.set_y(v); t.set_text(f"{v:.3e}")
        ax_cost_w.relim(); ax_cost_w.autoscale_view()

        # coûts non pondérés
        for b, v, t in zip(bars_u, C['u'], txt_u):
            b.set_height(v)
            t.set_y(v); t.set_text(f"{v:.3e}")
        ax_cost_u.relim(); ax_cost_u.autoscale_view()

        fig.canvas.draw_idle()

    for s in sliders:
        s.on_changed(update)

    def reset(_):
        for s in sliders:
            s.reset()

    btn_reset.on_clicked(reset)

    plt.show()
    return fig
