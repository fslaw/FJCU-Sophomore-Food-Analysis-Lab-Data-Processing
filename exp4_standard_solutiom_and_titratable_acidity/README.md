# 滴定酸度分析工具 (Titration Analysis Tools)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

這個專案提供了一套完整的滴定曲線分析工具，特別適用於食品和飲料樣品的滴定酸度分析。專案包含標準滴定曲線分析和進階數值微分分析兩個主要模組。

## 功能特點

### 1. 標準滴定曲線分析 (`titration_curve_en_v2.0.py`)
- 自動化數據處理和分析流程
- 精確計算滴定酸度
- 生成專業的滴定曲線圖表
- 提供詳細的實驗數據驗證
- 自動標注關鍵數據點

### 2. 進階微分分析 (`cubic spline微分模型.py`)
- 使用三次樣條插值進行曲線擬合
- 計算一階和二階導數
- 自動識別轉折點和 pKa 值
- 精確定位 pH 8.2 的滴定體積
- 輸出分段函數係數表格

## 系統需求

- Python 3.8+
- 相依套件：
  - numpy >= 1.25.2
  - pandas >= 2.1.0
  - matplotlib >= 3.7.2
  - scipy >= 1.13.1

## 快速開始

1. 複製專案：
```bash
git clone https://github.com/yourusername/titration-analysis.git
cd titration-analysis
```

2. 安裝相依套件：
```bash
pip install -r requirements.txt
```

3. 準備數據檔案：
- 將實驗數據整理為指定格式的 CSV 檔案
- 檔名：`standardization_and_titratable_acidity2_en.csv`
- 將檔案放在專案根目錄

4. 執行分析：
```bash
# 執行標準滴定曲線分析
python titration_curve_en_v2.0.py

# 執行進階微分分析
python "cubic spline微分模型.py"
```

## 輸入檔案格式

CSV 檔案應包含以下區塊：
1. KHP 標定區塊（標頭從第 4 行開始）
2. HCl 標定區塊（標頭從第 12 行開始）
3. 滴定曲線數據（標頭從第 34 行開始）

詳細格式請參考範例檔案。

## 輸出檔案

### 標準分析模式
- `titration_curve_plot_with_points.png`：滴定曲線圖
  - 包含誤差棒
  - 標注終點和半當量點
  - 顯示滴定酸度計算結果

### 進階分析模式
- `titration_curve_with_model_table.png`：進階分析圖
  - 包含一階和二階導數曲線
  - 標注轉折點和 pKa 值
  - 顯示 pH 8.2 位置
- `spline_piecewise_table.csv`：分段函數係數表

## 分析結果說明

### 1. 標準分析
- 計算並驗證 NaOH 和 HCl 的標準化結果
- 使用線性插值找到 pH 8.2 的滴定體積
- 計算滴定酸度（以檸檬酸計）

### 2. 進階分析
- 使用三次樣條提供更精確的曲線擬合
- 通過導數分析找到：
  - 當量點（二階導數為零）
  - pKa 值（一階導數極大值）
  - pH 8.2 精確位置

## 程式碼架構

```
titration-analysis/
├── titration_curve_en_v2.0.py    # 標準分析主程式
├── cubic spline微分模型.py        # 進階分析主程式
├── requirements.txt              # 相依套件列表
├── README.md                     # 說明文件
└── example_data/                 # 範例數據
    └── standardization_and_titratable_acidity2_en.csv
```

## 注意事項

- 確保 CSV 檔案格式正確
- 留意數據點的品質和完整性
- 注意實驗誤差的標示和計算

## 作者
T. Ong

## 授權

本專案採用 MIT 授權 - 詳見 [LICENSE](LICENSE) 檔案
