import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib

# Ensure English font
matplotlib.rcParams['font.family'] = 'DejaVu Sans'
matplotlib.rcParams['font.sans-serif'] = ['DejaVu Sans']

# --- File settings ---
CSV_FILE_NAME = 'standardization_and_titratable_acidity2_en.csv'

# --- Experimental constants ---
KHP_MW = 204.288  # g/mol
ACID_EQ_WT = 64.04  # equivalent weight (assumed citric acid)
SAMPLE_VOLUME_ML = 20.0  # sample volume in mL

# --- Settings specific to this script ---
ANNOTATE_POINTS = True  # Annotate each data point with (x, y)
POINT_ANNOTATION_FMT = "({x:.3f}, {y:.2f})"  # format for point labels
# Reduce font size to half for annotations (previously 8 -> now 4)
POINT_ANNOTATION_FONTSIZE = 4


def load_data(file_name):
    try:
        df_naoh = pd.read_csv(
            file_name,
            header=3,
            nrows=2,
            usecols=[0, 1, 4, 5],
            encoding='utf-8'
        )
        df_naoh.columns = ['Replicate', 'KHP mass (g)', 'Volume NaOH titrated (mL)', 'N NaOH']

        df_hcl = pd.read_csv(
            file_name,
            header=11,
            nrows=2,
            usecols=[0, 1, 2, 3],
            encoding='utf-8'
        )
        df_hcl.columns = ['Replicate', 'HCl volume (mL)', 'Volume NaOH titrated (mL)', 'N HCl']

        df_titration = pd.read_csv(
            file_name,
            header=33,
            usecols=[0, 1],
            encoding='utf-8'
        )
        df_titration.columns = ['mL NaOH', 'pH']
        df_titration = df_titration.dropna().astype(float)

        print(f"Successfully read all experimental data from '{file_name}'.")
        return df_naoh, df_hcl, df_titration
    except FileNotFoundError:
        print(f"Error: file '{file_name}' not found.")
        print("Please place the CSV in the same folder as this script.")
        return None, None, None
    except Exception as e:
        print(f"Error reading CSV: {e}")
        print("HINT: Check if the CSV structure matches expectations.")
        return None, None, None


def verify_calculations(df_naoh, df_hcl):
    print("\n--- A.1: Verify standardization of 0.1N NaOH ---")
    khp_col = 'KHP mass (g)'
    vol_col = 'Volume NaOH titrated (mL)'
    rep_col = 'Replicate'
    n_naoh_col = 'N NaOH'

    df_naoh['N NaOH (calculated)'] = (df_naoh[khp_col] / KHP_MW) / (df_naoh[vol_col] / 1000)
    print(df_naoh[[rep_col, n_naoh_col, 'N NaOH (calculated)']].to_string())

    calc_mean_naoh = df_naoh['N NaOH (calculated)'].mean()
    calc_std_naoh = df_naoh['N NaOH (calculated)'].std()
    csv_mean_naoh = 0.09910
    csv_std_naoh = 0.0002794

    print(f"\nYour calculated N NaOH Mean (from CSV): {csv_mean_naoh}")
    print(f"Script verification N NaOH Mean (x̄):   {calc_mean_naoh:.5f}")
    print(f"\nYour calculated N NaOH Stdev (from CSV): {csv_std_naoh}")
    print(f"Script verification N NaOH Stdev (SD):   {calc_std_naoh:.7f}")

    print("\n--- A.2: Verify standardization of 0.1N HCl ---")
    N_NaOH_MEAN_FROM_CSV = csv_mean_naoh
    rep_hcl = 'Replicate'
    v_naoh_col = 'Volume NaOH titrated (mL)'
    v_hcl_col = 'HCl volume (mL)'
    n_hcl_col = 'N HCl'

    df_hcl['N HCl (calculated)'] = (N_NaOH_MEAN_FROM_CSV * df_hcl[v_naoh_col]) / df_hcl[v_hcl_col]
    print(df_hcl[[rep_hcl, n_hcl_col, 'N HCl (calculated)']].to_string())

    calc_mean_hcl = df_hcl['N HCl (calculated)'].mean()
    calc_std_hcl = df_hcl['N HCl (calculated)'].std()
    csv_mean_hcl = 0.09910
    csv_std_hcl = 0.003504

    print(f"\nYour calculated N HCl Mean (from CSV): {csv_mean_hcl}")
    print(f"Script verification N HCl Mean (x̄):   {calc_mean_hcl:.5f}")
    print(f"\nYour calculated N HCl Stdev (from CSV): {csv_std_hcl}")
    print(f"Script verification N HCl Stdev (SD):   {calc_std_hcl:.7f}")

    return N_NaOH_MEAN_FROM_CSV


def analyze_titration_curve(df_titration, N_naoh):
    print("\n--- B.2: Analyze titration curve (Sample C) ---")

    x_volume = df_titration['mL NaOH']
    y_ph = df_titration['pH']

    try:
        V_endpoint = np.interp(8.2, y_ph, x_volume)
        print(f"Endpoint volume at pH 8.2 (interpolated): {V_endpoint:.3f} mL")
    except Exception as e:
        V_endpoint = x_volume[y_ph.sub(8.2).abs().idxmin()]
        print(f"Could not compute volume at pH 8.2 (error: {e}). Using nearest point instead: V = {V_endpoint} mL")

    numerator = V_endpoint * N_naoh * ACID_EQ_WT
    denominator = SAMPLE_VOLUME_ML * 10
    percent_acid = numerator / denominator
    print(f"Titratable acidity (as citric acid equivalent): {percent_acid:.4f} %")

    V_half_endpoint = V_endpoint / 2.0
    pKa_estimated = np.interp(V_half_endpoint, x_volume, y_ph)
    print(f"Half-equivalence point volume: {V_half_endpoint:.3f} mL")
    print(f"Estimated major pKa (pH at half-equivalence): {pKa_estimated:.2f}")

    propagated_absolute_error = 0.068
    print(f"Adding fixed absolute error (from propagation): +/- {propagated_absolute_error} mL to all points.")

    plt.figure(figsize=(10, 7))

    plt.errorbar(
        x_volume,
        y_ph,
        xerr=propagated_absolute_error,
        marker='o',
        linestyle='-',
        label=f'Titration data (error ± {propagated_absolute_error:.3f} mL)',
        capsize=5
    )

    # Mark endpoint and half-equivalence
    plt.axhline(y=8.2, color='red', linestyle='--', label='Endpoint (pH 8.2)')
    plt.axvline(x=V_endpoint, color='red', linestyle='--', label=f'Endpoint volume ({V_endpoint:.3f} mL)')
    plt.plot(V_endpoint, 8.2, 'r*', markersize=15)

    plt.axhline(y=pKa_estimated, color='green', linestyle=':', label=f'Estimated pKa ({pKa_estimated:.2f})')
    plt.axvline(x=V_half_endpoint, color='green', linestyle=':', label=f'Half-equivalence volume ({V_half_endpoint:.3f} mL)')
    plt.plot(V_half_endpoint, pKa_estimated, 'go', markersize=10)

    # Annotate each data point with its (x, y)
    if ANNOTATE_POINTS:
        # Smart placement: try a list of candidate offsets (in offset-points) and pick the first
        # location that does not overlap previously placed annotation boxes.
        fig = plt.gcf()
        ax = plt.gca()
        canvas = fig.canvas
        # ensure renderer exists
        canvas.draw()
        renderer = canvas.get_renderer()

        placed_bboxes = []  # list of Bbox in display coords
        xv_list = list(x_volume)
        yv_list = list(y_ph)
        n = len(xv_list)

        # candidate offsets (dx, dy, ha, va)
        candidates = [
            (6, 6, 'left', 'bottom'),
            (6, -6, 'left', 'top'),
            (-6, 6, 'right', 'bottom'),
            (-6, -6, 'right', 'top'),
            (12, 0, 'left', 'center'),
            (-12, 0, 'right', 'center'),
            (0, 12, 'center', 'bottom'),
            (0, -12, 'center', 'top')
        ]

        for i, (xv, yv) in enumerate(zip(xv_list, yv_list)):
            label = POINT_ANNOTATION_FMT.format(x=xv, y=yv)
            placed = False

            # transform data point to display coords
            disp_pt = ax.transData.transform((xv, yv))

            for dx, dy, ha, va in candidates:
                ann = plt.annotate(
                    label,
                    xy=(xv, yv),
                    xytext=(dx, dy),
                    textcoords='offset points',
                    fontsize=POINT_ANNOTATION_FONTSIZE,
                    ha=ha,
                    va=va,
                    bbox=dict(facecolor='white', alpha=0.9, edgecolor='black', linewidth=0.3, boxstyle='round,pad=0.2'),
                    zorder=5
                )

                # draw and get bbox in display coords
                try:
                    canvas.draw()
                    bb = ann.get_window_extent(renderer=renderer)
                except Exception:
                    # if renderer not available, accept default
                    bb = None

                conflict = False
                if bb is not None:
                    # check overlap with previously placed bboxes
                    for existing in placed_bboxes:
                        # small padding
                        if bb.overlaps(existing):
                            conflict = True
                            break

                    # also avoid the label covering the marker itself (ensure bbox does not contain the data point)
                    if not conflict:
                        x_px, y_px = disp_pt
                        if bb.contains(x_px, y_px):
                            conflict = True

                if conflict:
                    ann.remove()
                    continue
                else:
                    # accept annotation
                    if bb is not None:
                        placed_bboxes.append(bb)
                    placed = True
                    break

            if not placed:
                # fallback: place above-right with small offset
                ann = plt.annotate(
                    label,
                    xy=(xv, yv),
                    xytext=(6, 6),
                    textcoords='offset points',
                    fontsize=POINT_ANNOTATION_FONTSIZE,
                    ha='left',
                    va='bottom',
                    bbox=dict(facecolor='white', alpha=0.9, edgecolor='black', linewidth=0.3, boxstyle='round,pad=0.2'),
                    zorder=5
                )
                try:
                    canvas.draw()
                    bb = ann.get_window_extent(renderer=renderer)
                    if bb is not None:
                        placed_bboxes.append(bb)
                except Exception:
                    pass

    plt.title('pH Titration Curve of Apple Juice (Sample C)', fontsize=16)
    plt.xlabel(f'Used volume (mL) of {N_naoh}N NaOH(aq)', fontsize=12)
    plt.ylabel('pH value', fontsize=12)
    # Place legend in upper-left but ensure it has a solid frame and is on top
    lg = plt.legend(loc='upper left', frameon=True, fontsize=9)
    if lg is not None:
        lg.get_frame().set_alpha(0.9)
        try:
            lg.set_zorder(20)
            lg.get_frame().set_linewidth(0.5)
        except Exception:
            pass
    plt.grid(True)

    plot_filename = 'titration_curve_plot_with_points.png'
    # on-plot annotation for fixed x-error
    # (Error textbox removed by user request)

    plt.savefig(plot_filename, dpi=300)
    print(f"\nPlot saved as '{plot_filename}' (300 dpi)")

    plt.show()


if __name__ == '__main__':
    if not os.path.exists(CSV_FILE_NAME):
        print(f"Error: could not find file '{CSV_FILE_NAME}'.")
        print(f"Please make sure '{CSV_FILE_NAME}' is in the script folder.")
    else:
        df_naoh, df_hcl, df_titration = load_data(CSV_FILE_NAME)
        if df_naoh is not None and df_hcl is not None and df_titration is not None:
            N_naoh_mean = verify_calculations(df_naoh, df_hcl)
            analyze_titration_curve(df_titration, N_naoh_mean)
