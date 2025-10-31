# 滴定酸度分析工具 (Titration Analysis Tools)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

這個專案提供了一套完整的滴定曲線分析工具，特別適用於食品和飲料樣品的滴定酸度分析。專案包含標準滴定曲線分析（SOP）和兩種進階數值微分分析模組。

## 功能特點

### 1. 標準滴定曲線分析 (`titration_curve_en_v2.0.py`)
- **SOP 遵循**：嚴格依據實驗手冊，使用**線性內插法** (`numpy.interp`) 找出 **$pH 8.2$** (酚酞終點) 所對應的「滴定終點」體積 (V = 9.343 mL)。
- **結果計算**：自動計算滴定酸度 (TA = 0.2965 %) 與表觀 $pKa$ (4.87)。
- **誤差模型**：包含基於誤差疊加原理的**固定絕對誤差線** (±0.068 mL)。

### 2. 進階微分分析 (Cubic Spline) (`cubic spline微分模型.py`)
- **模型建立**：使用**三次樣條插值** (`scipy.interpolate.CubicSpline`) 對**原始數據** 建立平滑函數。
- **問題展示**：計算一階和二階導數，展示因**過度擬合 (Overfitting)** 實驗雜訊 而產生的 **7+ 個「偽特徵」拐點**。
- **模型匯出**：輸出構成曲線的 10 條分段函數係數至 `spline_piecewise_table.csv`。

### 3. 進階微分分析 (Savitzky-Golay 濾波) (`濾波平滑微分.py`)
- **訊號處理**：(推薦) 此腳本 首先使用 **Savitzky-Golay 濾波器** 對原始 $pH$ 數據 進行**平滑化 (Smoothing)**，以去除「短週期雜訊」。
- **精確分析**：接著，對**平滑後的數據** 進行 Cubic Spline 與二階導數分析。
- **真實特徵辨識**：成功濾除雜訊，辨識出**化學特徵點**：
    - **真實當量點** (d²=0, d' max): V = **8.502 mL** (at pH 7.02)
    - **真實 pKa** (d²=0, d' min): V = **2.508 mL** (at pH 4.41)
- **結果比較**：清楚地揭示了「滴定終點」(9.34 mL) 與「化學當量點」(8.50 mL) 之間的顯著差異。

## 系統需求

- Python 3.8+
- 相依套件 (詳見 `requirements.txt`)：
  - numpy
  - pandas
  - matplotlib
  - scipy

## 快速開始

1. 複製專案：
```bash
git clone [https://github.com/yourusername/FJCU-Food-Science-Sophomore-Lab-Data.git](https://github.com/yourusername/FJCU-Food-Science-Sophomore-Lab-Data.git)
cd FJCU-Food-Science-Sophomore-Lab-Data/Exp4_Standard_Solutions_and_Titratable_Acidity

## 授權

本專案採用 MIT 授權 - 詳見 [LICENSE](LICENSE) 檔案
