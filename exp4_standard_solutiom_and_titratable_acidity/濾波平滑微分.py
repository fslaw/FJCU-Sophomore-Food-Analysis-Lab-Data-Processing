import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.signal import savgol_filter

# --- 檔案設定 ---
CSV_FILE_NAME = 'standardization_and_titratable_acidity2_en.csv'

# --- 實驗參數 ---
N_naoh = 0.0991  # 經標定後的 NaOH 溶液濃度

# --- 平滑化參數 ---
SAVGOL_WINDOW = 5 
SAVGOL_POLY_ORDER = 2

# --- 1. 載入滴定數據 ---
try:
    df_titration = pd.read_csv(
        CSV_FILE_NAME,
        header=33,      # 標頭在第 34 行 (index 33)
        usecols=['mL NaOH', 'pH'], # 只讀取第一組 "mL NaOH" 和 "pH"
        encoding='utf-8'
    )
    df_titration = df_titration.dropna().astype(float)
    
    v_data = df_titration['mL NaOH'].values
    ph_data = df_titration['pH'].values
    
    print(f"成功從 '{CSV_FILE_NAME}' 載入 {len(v_data)} 個滴定數據點。")

except Exception as e:
    print(f"讀取 CSV 時發生錯誤: {e}")
    exit()

# --- 2. 數據平滑化 (Filtering) ---
try:
    ph_smooth = savgol_filter(ph_data, SAVGOL_WINDOW, SAVGOL_POLY_ORDER)
    print(f"已應用 Savitzky-Golay 濾波器 (Window={SAVGOL_WINDOW}, Order={SAVGOL_POLY_ORDER})。")
except ValueError as e:
    print(f"SavGol 濾波器參數錯誤: {e}")
    exit()

# --- 3. 建立平滑函數 (Cubic Spline) ---
spline = CubicSpline(v_data, ph_smooth)

# --- 4. 計算平滑後的導數 ---
spline_d1 = spline.derivative(1) # 一階導數
spline_d2 = spline.derivative(2) # 二階導數

# --- 5. (重構) 找出所有關鍵點 ---

def find_volume_at_ph(spline_func, target_ph, v_min, v_max, tolerance=1e-6):
    """
    使用二分法找到特定 pH 值對應的體積
    """
    def objective(v):
        return float(spline_func(v)) - target_ph
    
    ph_min = float(spline_func(v_min))
    ph_max = float(spline_func(v_max))
    
    if (target_ph - ph_min) * (target_ph - ph_max) > 0:
        print(f"錯誤: 目標 pH {target_ph} 不在數據範圍 [{ph_min:.2f}, {ph_max:.2f}] 內。")
        return None
    
    low_v, high_v = (v_min, v_max) if ph_min < ph_max else (v_max, v_min)
    
    while (high_v - low_v) > tolerance:
        mid_v = (low_v + high_v) / 2
        mid_ph = float(spline_func(mid_v))
        
        if abs(mid_ph - target_ph) < tolerance:
            return mid_v
        elif mid_ph < target_ph:
            low_v = mid_v
        else:
            high_v = mid_v
            
    return (low_v + high_v) / 2

print("\n--- (平滑後) 分析結果 ---")

# (1) 找出 pH 8.2 的位置
v_ph82 = find_volume_at_ph(spline, 8.2, v_data.min(), v_data.max())
if v_ph82 is not None:
    print(f"法定終點 (pH 8.2) 對應的體積: {v_ph82:.3f} mL")

# (2) 找出所有轉折點 (二階導數=0)
inflection_points_v = spline_d2.roots()
valid_inflection_points = []
for v in inflection_points_v:
    if v_data.min() < v < v_data.max():
        valid_inflection_points.append(v)
print(f"在數據範圍內找到 {len(valid_inflection_points)} 個轉折點 (d²(pH)/dV² = 0): {[f'{v:.3f}' for v in valid_inflection_points]} mL")

# (3) 找出 當量點 (EQ) (斜率最大)
v_eq = None
max_d1 = float('-inf')

for v in valid_inflection_points:
    d1 = float(spline_d1(v))
    if d1 > max_d1:
        max_d1 = d1
        v_eq = v

if v_eq is not None:
    ph_eq = float(spline(v_eq))
    print(f"真實當量點 (最大斜率) 位於: V = {v_eq:.3f} mL, pH = {ph_eq:.2f}")
else:
    print("未找到主要的當量點。")


# --- 6. 繪製圖表 (已修正標示) ---
v_smooth_axis = np.linspace(v_data.min(), v_data.max(), 500)
fig, ax1 = plt.subplots(figsize=(12, 8))

# a) 繪製 "原始" 數據點
ax1.plot(v_data, ph_data, 'ko', label='Original Data (Noisy)', markersize=8, zorder=10)

# b) 繪製 "平滑後" 的曲線
ax1.plot(v_smooth_axis, spline(v_smooth_axis), 'b-', label='Smoothed Curve (SavGol + Spline)', linewidth=3)

# --- 建立第二個 Y 軸 (用於導數) ---
ax2 = ax1.twinx()
ax2.plot(v_smooth_axis, spline_d1(v_smooth_axis), 'g--', label='1st Derivative (Smoothed)')
ax2.plot(v_smooth_axis, spline_d2(v_smooth_axis), 'r:', label='2nd Derivative (Smoothed)')
ax2.axhline(y=0, color='black', linestyle='-', linewidth=0.7)

# --- f) 標示所有關鍵點 ---

# (1) 標示 pH 8.2 的位置
if v_ph82 is not None:
    ax1.plot(v_ph82, 8.2, 'k^', markersize=12, zorder=12,
             label=f'Legal Endpoint (pH 8.2)\nV = {v_ph82:.3f} mL')
    ax1.axvline(x=v_ph82, color='k', linestyle='-.', linewidth=1)
    ax1.axhline(y=8.2, color='k', linestyle='-.', linewidth=1)

# (2) 標示所有轉折點 (pKa 或 EQ)
for v_ip in valid_inflection_points:
    ph_ip = float(spline(v_ip)) # 取得該點的 pH
    
    # 檢查這個點是否為我們找到的 "主要當量點"
    if v_eq is not None and abs(v_ip - v_eq) < 1e-3:
        # 如果是，標示為 Equivalence Point
        label = f'Equivalence Point (Inflection)\nV = {v_ip:.3f} mL, pH = {ph_ip:.2f}'
        ax1.plot(v_ip, ph_ip, 'r*', markersize=20, zorder=11, label=label)
        ax1.axvline(x=v_ip, color='r', linestyle=':', linewidth=2)
    else:
        # *** 依照您的要求 ***
        # 如果不是主要當量點，則一律標示為 pKa (Inflection)
        label = f'pKa (Inflection)\nV = {v_ip:.3f} mL, pH = {ph_ip:.2f}'
        ax1.plot(v_ip, ph_ip, 'mP', markersize=15, zorder=11, label=label) # P = Plus
        ax1.axvline(x=v_ip, color='m', linestyle=':', linewidth=2)

# g) 圖表美化
ax1.set_title('Derivative Analysis of Titration Curve (Smoothed)', fontsize=16)
ax1.set_xlabel(f'Used volume (mL) of {N_naoh:.4f}N NaOH(aq)', fontsize=12)
ax1.set_ylabel('pH value', fontsize=12)
ax2.set_ylabel('Derivative value', fontsize=12)

# 合併兩個 Y 軸的圖例
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper left', frameon=True, fontsize=9)
ax1.grid(True)
ax2.grid(False) # 關閉第二個 Y 軸的格線

plot_filename = 'titration_curve_savgol_analysis_v3.png'
plt.savefig(plot_filename, dpi=300)
print(f"\n平滑化導數分析圖表已儲存為 '{plot_filename}'")
plt.show()
