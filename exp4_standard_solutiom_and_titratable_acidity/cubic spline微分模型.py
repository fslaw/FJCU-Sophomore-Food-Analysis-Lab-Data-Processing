import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import CubicSpline  # 匯入三次樣條插值

# --- 1. 檔案設定 ---
CSV_FILE_NAME = 'standardization_and_titratable_acidity2_en.csv'

# 定義 NaOH 的實際標定濃度
N_naoh = 0.0991  # 經標定後的 NaOH 溶液濃度

# --- 2. 載入滴定數據 ---
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

def find_volume_at_ph(spline_func, target_ph, v_min, v_max, tolerance=1e-6):
    """
    使用二分法找到特定 pH 值對應的體積
    """
    def objective(v):
        return float(spline_func(v)) - target_ph
    
    # 檢查目標 pH 是否在搜尋範圍的 pH 內
    ph_min = float(spline_func(v_min))
    ph_max = float(spline_func(v_max))
    
    if (target_ph - ph_min) * (target_ph - ph_max) > 0:
        # 如果 target_ph 同時大於或小於範圍兩端，代表不在範圍內
        print(f"錯誤: 目標 pH {target_ph} 不在數據範圍 [{ph_min:.2f}, {ph_max:.2f}] 內。")
        return None
    
    # 二分法搜尋
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

# --- 3. 建立平滑函數 (Cubic Spline) ---
spline = CubicSpline(v_data, ph_data)

# --- 4. (新功能) 準備表格數據 ---
print("\n--- 準備 Cubic Spline (三次樣條) 區間函數表格 ---")
# 函數形式為: pH(V) = C3(V-V_i)^3 + C2(V-V_i)^2 + C1(V-V_i) + C0
table_data = []
table_cols = ['Interval (V_i to V_i+1)', 'C3 (x-V_i)^3', 'C2 (x-V_i)^2', 'C1 (x-V_i)', 'C0 (pH at V_i)']

for i in range(len(v_data) - 1):
    v_start = v_data[i]
    v_end = v_data[i+1]
    
    c3 = spline.c[0, i]
    c2 = spline.c[1, i]
    c1 = spline.c[2, i]
    c0 = spline.c[3, i] # 這等於 ph_data[i]
    
    # 將數據格式化為字串
    row = [
        f"[{v_start}, {v_end}]",
        f"{c3:.3e}",
        f"{c2:.3e}",
        f"{c1:.3e}",
        f"{c0:.2f}"
    ]
    table_data.append(row)
print("表格數據準備完畢。")
# ------------------------------------

# --- 5. 計算一階與二階導數 ---
spline_d1 = spline.derivative(1)
spline_d2 = spline.derivative(2)

# --- 6. 找出轉折點 (Inflection Points) ---
inflection_points_v = spline_d2.roots()

print("\n--- 導數分析結果 ---")
# 尋找 pH 8.2 的體積
v_ph82 = find_volume_at_ph(spline, 8.2, v_data.min(), v_data.max())
if v_ph82 is not None:
    print(f"pH 8.2 對應的體積: {v_ph82:.3f} mL")
else:
    print("無法找到 pH 8.2 對應的體積")

pKa_points_v = spline_d1.roots()
print(f"找到 pKa (緩衝點) 候選 位於: {pKa_points_v} mL")
print(f"找到轉折點 (二階導數=0) 位於: {inflection_points_v} mL")

valid_inflection_points = []
for v in inflection_points_v:
    if v > v_data.min() and v < v_data.max():
        valid_inflection_points.append(v)
print(f"在數據範圍內的轉折點: {valid_inflection_points} mL")


# --- 7. 繪製圖表 (複合圖) ---
fig, (ax1, ax2) = plt.subplots(
    2, 1, 
    figsize=(12, 16), # 增加總高度以容納表格
    gridspec_kw={'height_ratios': [3, 1]}
)

# --- (A) 繪製上方的滴定曲線圖 (ax1) ---
v_smooth = np.linspace(v_data.min(), v_data.max(), 500)

ax1.plot(v_data, ph_data, 'ko', label='Original Data Points', markersize=8, zorder=10)
ax1.plot(v_smooth, spline(v_smooth), 'b-', label='Titration Curve (Cubic Spline)', linewidth=2)
ax1.plot(v_smooth, spline_d1(v_smooth), 'g--', label='1st Derivative (d(pH)/dV)')
ax1.plot(v_smooth, spline_d2(v_smooth), 'r:', label='2nd Derivative (d²(pH)/dV²)')
ax1.axhline(y=0, color='black', linestyle='-', linewidth=0.7)

# f) 標示 pH 8.2 的位置
if v_ph82 is not None:
    ph_82 = float(spline(v_ph82)) # 幾乎會等於 8.2
    ax1.plot(v_ph82, ph_82, 'k^', markersize=15, zorder=11, 
             label=f'pH 8.2 Point\nV = {v_ph82:.3f} mL')
    ax1.axvline(x=v_ph82, color='k', linestyle='-.', linewidth=1)

# g) 標示轉折點 (*** 已修改 ***)
for v_ip in valid_inflection_points:
    ph_ip = float(spline(v_ip)) # 取得該點的 pH
    d1_ip = float(spline_d1(v_ip))
    
    # 建立標籤，包含 pH (Y值)
    label_text = f'V = {v_ip:.3f} mL, pH = {ph_ip:.2f}'
    
    if d1_ip > 1.0: # 當量點
        label = f'Equivalence Point (Inflection)\n{label_text}'
        ax1.plot(v_ip, ph_ip, 'r*', markersize=20, zorder=11, label=label)
        ax1.axvline(x=v_ip, color='r', linestyle=':', linewidth=2)
    else: # pKa
        label = f'pKa (Inflection)\n{label_text}'
        ax1.plot(v_ip, ph_ip, 'mP', markersize=15, zorder=11, label=label)
        ax1.axvline(x=v_ip, color='m', linestyle=':', linewidth=2)

# h) 圖表美化
ax1.set_title('pH Titration Curve of Apple Juice (Sample C)', fontsize=16)
ax1.set_xlabel(f'Used volume (mL) of {N_naoh:.4f}N NaOH(aq)', fontsize=12)
ax1.set_ylabel('pH value / Derivative value', fontsize=12)
ax1.legend(loc='upper left', frameon=True, fontsize=9)
ax1.grid(True)
ax1.set_ylim(min(ph_data.min(), spline_d2(v_smooth).min()) - 1, ph_data.max() + 1)

# --- (B) 繪製下方的表格 (ax2) ---
ax2.axis('off') # 關閉 ax2 的 x, y 軸
ax2.set_title("Cubic Spline Model Coefficients (Piecewise Functions)", pad=20, fontsize=14)
# 函數形式為: pH(V) = C3(V-V_i)^3 + C2(V-V_i)^2 + C1(V-V_i) + C0
tbl = ax2.table(
    cellText=table_data,
    colLabels=table_cols,
    loc='center',
    cellLoc='center'
)
tbl.auto_set_font_size(False)
tbl.set_fontsize(9)
tbl.scale(1.0, 1.4) # 調整表格大小 (寬度, 高度)

# --- 8. 儲存並顯示 ---
plot_filename = 'titration_curve_with_model_table.png'
plt.tight_layout() # 自動調整排版
plt.savefig(plot_filename, dpi=300)
print(f"\n導數分析圖表 (含模型) 已儲存為 '{plot_filename}'")
plt.show()

# 儲存表格數據到 CSV 檔案
def save_table_to_csv():
    """將區間函數表格數據儲存為 CSV 檔案"""
    df = pd.DataFrame(table_data, columns=table_cols)
    csv_filename = 'spline_piecewise_table.csv'
    df.to_csv(csv_filename, index=False)
    print(f"\n區間函數表格已儲存為 '{csv_filename}'")

# --- 9. 主執行區 ---
if __name__ == '__main__':
    # 檢查必要的檔案是否存在
    if not os.path.exists(CSV_FILE_NAME):
        print(f"錯誤：找不到檔案 '{CSV_FILE_NAME}'")
        print(f"請確認 '{CSV_FILE_NAME}' 在腳本所在的資料夾中。")
        exit(1)
        
    try:
        print("開始執行滴定曲線分析...")
        # 資料載入與分析在腳本主體中已完成
        # 額外儲存表格數據
        save_table_to_csv()
        print("\n分析完成！")
    except Exception as e:
        print(f"執行過程中發生錯誤：{e}")
        exit(1)
        
