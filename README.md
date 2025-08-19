
This repository contains a Python pipeline for analyzing **time–density curves (TDCs)** in intracranial aneurysm imaging.  
It combines **SVD-based deconvolution** with **gamma variate fitting** to estimate perfusion parameters, impulse response functions (IRFs), and evaluate the effect of patient motion or noise.

---

## 🔑 Features
- **Time–Density Curve (TDC) extraction** from `.raw` angiographic projection data using aneurysm and inlet ROI masks (`.tif`).
- **API_parameters**:
  - Bolus Arrival Time (BAT)
  - Peak Height (PH)
  - Time-to-Peak (TTP)
  - Area Under Curve (AUC)
  - Mean Transit Time (MTT)
  - Max derivative
- **API_SVD**: Toeplitz matrix–based SVD deconvolution with Tikhonov regularization.
- **Gamma variate modeling**:
  - Fits gamma curves to inlet and aneurysm TDCs
  - Produces smoother estimates for parameter extraction
- **Patient motion simulation**: Demonstrates impact of noise and drift on TDCs.
- **Excel integration**: Writes computed parameters to specific columns in `.xlsx` files using `openpyxl`.
- **Visualization**:
  - TDCs (raw, noisy, gamma-fitted)
  - Impulse Response Function (IRF)
  - Convolved inlet vs aneurysm TDC
  - ROI mask overlays on projection images

---
