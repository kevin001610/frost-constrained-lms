# frost-constrained-lms

原始資料來源：O. L. Frost III, "An Algorithm for Linearly Constrained Adaptive Array Processing," *Proc. IEEE*, vol. 60, no. 8, pp. 926–935, Aug. 1972（張順雄教授 聲納訊號處理課程 論文報告）。

互動網頁：https://kevin001610.github.io/frost-constrained-lms/

- `index.html`：用 JavaScript 重新實作 MATLAB 復現程式的互動圖解（寬頻陣列、R_xx 重建、最佳受限權重、受限 LMS 與梯度投影法、輸出功率依來源拆解與空間-頻率響應、平均權重收斂、失調量與約束誤差的集成驗證）。
- `matlab/`：MATLAB 復現程式。先執行 `frost1972_verify.m`（24 項自我驗證），再執行 `frost1972_run.m` 產生 `matlab/results/` 的圖檔與數值摘要。

注意：`frost1972_adapt.m` 的 `quant()` 用 MATLAB `round`（遇到剛好一半時遠離零進位），6 位數算術裡這種「剛好一半」很常見，會讓梯度投影法的誤差出現系統性漂移；網頁的集成驗證可以切換捨入規則比較。
