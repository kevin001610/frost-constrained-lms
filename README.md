# frost-constrained-lms

原始資料來源：O. L. Frost III, "An Algorithm for Linearly Constrained Adaptive Array Processing," *Proc. IEEE*, vol. 60, no. 8, pp. 926–935, Aug. 1972（張順雄教授 聲納訊號處理課程 論文報告）。

互動網頁：https://kevin001610.github.io/frost-constrained-lms/

- `index.html`：用 JavaScript 重新實作 MATLAB 復現程式，只重現論文 Sec. VII 做過的模擬：Table I、Fig. 9（看向方向響應）、Fig. 10（入射頻譜）、R_xx 特徵值範圍、式 (26) 的 μ 上界與式 (27) 的失調量，以及 Fig. 11、Fig. 12（6 位有效數字、μ = 0.01、200 次更新的約束 LMS 與梯度投影法）。預設值就是論文的設定；Table I 的源參數、𝓕、μ、更新次數 N 與算術精度都可以在圖旁自行調整，並有按鈕回到論文設定。
- `matlab/`：MATLAB 復現程式。先執行 `frost1972_verify.m`（24 項自我驗證），再執行 `frost1972_run.m` 產生 `matlab/results/` 的圖檔與數值摘要。

