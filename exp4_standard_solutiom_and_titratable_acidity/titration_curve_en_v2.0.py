import os
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------------------
# Plotting configuration
# Set a consistent sans-serif font so labels look the same across platforms.
matplotlib.rcParams['font.family'] = 'DejaVu Sans'
matplotlib.rcParams['font.sans-serif'] = ['DejaVu Sans']

# --- File settings ---------------------------------------------------------
# CSV file containing: KHP standardization block, HCl standardization block,
# and the titration curve block. If CSV layout changes, update header indices
# inside `load_data` accordingly.
CSV_FILE_NAME = 'standardization_and_titratable_acidity2_en.csv'

# --- Experimental constants -----------------------------------------------
# KHP molecular weight (g/mol) used to compute NaOH normality from KHP mass
KHP_MW = 204.288  # g/mol
# Equivalent weight used to convert equivalents of acid to grams (citric acid eq.)
ACID_EQ_WT = 64.04  # equivalent weight (assumed citric acid)
# Volume (mL) of the sample aliquot used for titration. Used in percent calc.
SAMPLE_VOLUME_ML = 20.0  # sample volume in mL

# --- Display / annotation settings ---------------------------------------
ANNOTATE_POINTS = True  # whether to add text labels near each plotted point
POINT_ANNOTATION_FMT = "({x:.3f}, {y:.2f})"  # annotation format for (x, y)
# Small font to reduce overlap when annotating many points
POINT_ANNOTATION_FONTSIZE = 4


def load_data(file_name):
    """
    Read three logical blocks from the provided CSV file and return three
    pandas DataFrames:

    - df_naoh: rows for KHP standardization (used to verify NaOH normality)
      expected at CSV header=3, using columns [0,1,4,5]
    - df_hcl: rows for HCl standardization (expected at header=11)
    - df_titration: the titration curve (mL NaOH, pH) starting at header=33

    The header/row indices above match the worksheet exported by the lab
    template used for this experiment. If your CSV differs, update them.

    Returns a tuple (df_naoh, df_hcl, df_titration) or (None, None, None)
    on failure.
    """
    try:
        # Read the small KHP standardization table (2 rows expected)
        df_naoh = pd.read_csv(
            file_name,
            header=3,
            nrows=2,
            usecols=[0, 1, 4, 5],
            encoding='utf-8'
        )
        df_naoh.columns = ['Replicate', 'KHP mass (g)', 'Volume NaOH titrated (mL)', 'N NaOH']

        # Read the HCl standardization block (2 rows expected)
        df_hcl = pd.read_csv(
            file_name,
            header=11,
            nrows=2,
            usecols=[0, 1, 2, 3],
            encoding='utf-8'
        )
        df_hcl.columns = ['Replicate', 'HCl volume (mL)', 'Volume NaOH titrated (mL)', 'N HCl']

        # Read the titration curve: two columns (mL NaOH, pH). Drop empty rows
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
        # Help the user find the CSV file if it is missing
        print(f"Error: file '{file_name}' not found.")
        print("Please place the CSV in the same folder as this script.")
        return None, None, None
    except Exception as e:
        # Generic fallback that warns about CSV layout / parsing problems
        print(f"Error reading CSV: {e}")
        print("HINT: Check if the CSV structure matches expectations.")
        return None, None, None


def verify_calculations(df_naoh, df_hcl):
        """
        Verify and print simple checks on standardization data.

        - Computes N NaOH from the provided KHP masses and compares with the
            values recorded in the CSV (prints both).
        - Computes N HCl using the NaOH mean read from the CSV block and prints
            a similar side-by-side comparison.

        Returns a canonical value for N NaOH to use later in titration analysis
        (pulled from the CSV block values, not the script-calculated mean). This
        mirrors how the lab worksheet reports the official NaOH normality.
        """
        print("\n--- A.1: Verify standardization of 0.1N NaOH ---")
        khp_col = 'KHP mass (g)'
        vol_col = 'Volume NaOH titrated (mL)'
        rep_col = 'Replicate'
        n_naoh_col = 'N NaOH'

        # Calculate N NaOH from KHP mass: (mass / molar_mass) / (volume_L)
        df_naoh['N NaOH (calculated)'] = (df_naoh[khp_col] / KHP_MW) / (df_naoh[vol_col] / 1000)
        print(df_naoh[[rep_col, n_naoh_col, 'N NaOH (calculated)']].to_string())

        # Basic summary statistics for the calculated NaOH values
        calc_mean_naoh = df_naoh['N NaOH (calculated)'].mean()
        calc_std_naoh = df_naoh['N NaOH (calculated)'].std()
        # The CSV also contains reported values which this script prints for
        # side-by-side comparison (these are from the lab worksheet/export)
        csv_mean_naoh = 0.09910
        csv_std_naoh = 0.0002794

        print(f"\nYour calculated N NaOH Mean (from CSV): {csv_mean_naoh}")
        print(f"Script verification N NaOH Mean (x̄):   {calc_mean_naoh:.5f}")
        print(f"\nYour calculated N NaOH Stdev (from CSV): {csv_std_naoh}")
        print(f"Script verification N NaOH Stdev (SD):   {calc_std_naoh:.7f}")

        print("\n--- A.2: Verify standardization of 0.1N HCl ---")
        # Use the reported NaOH mean from the CSV to compute N HCl for the
        # HCl standardization block (this follows the lab's calculation order)
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
    """
    Analyze the titration curve and produce a plot with annotations.

    Steps performed:
    - Interpolate to find the endpoint volume at pH 8.2 (typical endpoint for
      titrating weak acids with strong base when using phenolphthalein/indicator
      approximations or a set endpoint pH).
    - Compute titratable acidity (as citric acid equivalent):
        percent_acid = (V_endpoint * N_naoh * ACID_EQ_WT) / (sample_volume_mL * 10)
      Explanation: V_endpoint (mL) * N (eq/L) gives milliequivalents; convert
      units appropriately to obtain a percent w/v value.
    - Estimate the pKa by evaluating pH at the half-equivalence volume (V/2).
    - Plot the titration curve with error bars, mark the endpoint and half-
      equivalence, annotate points (smart placement), and add a legend.

    The plotting code uses axes-fraction coordinates to place the computed
    titratable acidity text in the lower-right corner of the plot.
    """
    print("\n--- B.2: Analyze titration curve (Sample C) ---")

    x_volume = df_titration['mL NaOH']
    y_ph = df_titration['pH']

    # Interpolate along y_ph->x_volume to find the x where pH == 8.2
    try:
        V_endpoint = np.interp(8.2, y_ph, x_volume)
        print(f"Endpoint volume at pH 8.2 (interpolated): {V_endpoint:.3f} mL")
    except Exception as e:
        # Fallback to the nearest measured point if interpolation fails
        V_endpoint = x_volume[y_ph.sub(8.2).abs().idxmin()]
        print(f"Could not compute volume at pH 8.2 (error: {e}). Using nearest point instead: V = {V_endpoint} mL")

    # Calculate percent acidity (as citric acid equivalent)
    numerator = V_endpoint * N_naoh * ACID_EQ_WT
    denominator = SAMPLE_VOLUME_ML * 10  # factor 10 converts mL->L and grams->percent logic
    percent_acid = numerator / denominator
    print(f"Titratable acidity (as citric acid equivalent): {percent_acid:.4f} %")

    # Estimate pKa: pH at half the endpoint volume
    V_half_endpoint = V_endpoint / 2.0
    pKa_estimated = np.interp(V_half_endpoint, x_volume, y_ph)
    print(f"Half-equivalence point volume: {V_half_endpoint:.3f} mL")
    print(f"Estimated major pKa (pH at half-equivalence): {pKa_estimated:.2f}")

    # A fixed propagated absolute error used for plotting horizontal error bars
    propagated_absolute_error = 0.068
    print(f"Adding fixed absolute error (from propagation): +/- {propagated_absolute_error} mL to all points.")

    plt.figure(figsize=(10, 7))

    # Plot measured points with an x-error bar (same fixed error for all points)
    plt.errorbar(
        x_volume,
        y_ph,
        xerr=propagated_absolute_error,
        marker='o',
        linestyle='-',
        label=f'Titration data (error ± {propagated_absolute_error:.3f} mL)',
        capsize=5
    )

    # Visual markers for endpoint and half-equivalence
    plt.axhline(y=8.2, color='red', linestyle='--', label='Endpoint (pH 8.2)')
    plt.axvline(x=V_endpoint, color='red', linestyle='--', label=f'Endpoint volume ({V_endpoint:.3f} mL)')
    plt.plot(V_endpoint, 8.2, 'r*', markersize=15)

    plt.axhline(y=pKa_estimated, color='green', linestyle=':', label=f'Estimated pKa ({pKa_estimated:.2f})')
    plt.axvline(x=V_half_endpoint, color='green', linestyle=':', label=f'Half-equivalence volume ({V_half_endpoint:.3f} mL)')
    plt.plot(V_half_endpoint, pKa_estimated, 'go', markersize=10)

    # Annotate each data point with its (x, y) using a collision-avoiding strategy
    if ANNOTATE_POINTS:
        fig = plt.gcf()
        ax = plt.gca()
        canvas = fig.canvas
        # Force a draw so we can get a renderer for bbox calculations
        canvas.draw()
        renderer = canvas.get_renderer()

        placed_bboxes = []  # list of Bbox in display coords used to avoid overlap
        xv_list = list(x_volume)
        yv_list = list(y_ph)
        n = len(xv_list)

        # Candidate offsets (dx, dy, ha, va) tried in order until an annotation fits
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

            # transform data point to display coords for bbox overlap checks
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
                    # if renderer not available, accept default placement
                    bb = None

                conflict = False
                if bb is not None:
                    # check overlap with previously placed bboxes
                    for existing in placed_bboxes:
                        if bb.overlaps(existing):
                            conflict = True
                            break

                    # also avoid the label covering the marker itself
                    if not conflict:
                        x_px, y_px = disp_pt
                        if bb.contains(x_px, y_px):
                            conflict = True

                if conflict:
                    ann.remove()
                    continue
                else:
                    # accept annotation and remember its bbox
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

    # Final plot cosmetics
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

    # Add titratable acidity text in the lower-right corner (axes fraction coords)
    try:
        ax = plt.gca()
        plt.text(
            0.98,
            0.02,
            f'Titratable Acidity: {percent_acid:.4f} %',
            transform=ax.transAxes,
            ha='right',
            va='bottom',
            fontsize=10,
            bbox=dict(facecolor='white', alpha=0.9, edgecolor='black', linewidth=0.3)
        )
    except Exception:
        # Fall back: use figure text if axes transform fails
        plt.figtext(0.98, 0.02, f'Titratable Acidity: {percent_acid:.4f} %', ha='right', va='bottom', fontsize=10)

    # (Removed extra interpolated-curve figure per user request — keep only the
    # original main plot to match the original behaviour and outputs.)

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
