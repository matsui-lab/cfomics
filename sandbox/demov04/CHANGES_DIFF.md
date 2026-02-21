# DR-Learner 修正内容（diff 形式）

## ファイル: `packages/cfomics/R/methods_drlearner.R`

### 変更1: 関数シグネチャに新規引数を追加

```diff
 cf_fit_drlearner <- function(X, T, Y,
                              covariate_names = NULL,
                              random_state = 0L,
                              model_propensity = "auto",
                              model_regression = "auto",
+                             model_final = "auto",
                              cv = 5L,
+                             min_propensity = 0.05,
+                             standardize_y = TRUE,
                              ...) {
```

---

### 変更2: Y の標準化処理を追加（FIX #4）

```diff
   econml <- reticulate::import("econml")
+  sklearn <- reticulate::import("sklearn")
   np <- reticulate::import("numpy")
   DRLearner <- reticulate::import("econml.dr")$DRLearner

+  # Y の標準化
+  Y_orig <- Y
+  if (standardize_y) {
+    Y_mean <- mean(Y)
+    Y_sd <- sd(Y)
+    Y <- (Y - Y_mean) / Y_sd
+  } else {
+    Y_mean <- 0
+    Y_sd <- 1
+  }
+
   if (!is.null(random_state)) {
     np$random$seed(as.integer(random_state))
+    random <- reticulate::import("random")
+    random$seed(as.integer(random_state))
   }
```

---

### 変更3: デフォルトモデルの設定（FIX #1）

```diff
+  # Propensity model
+  if (identical(model_propensity, "auto")) {
+    LogisticRegressionCV <- sklearn$linear_model$LogisticRegressionCV
+    model_prop_py <- LogisticRegressionCV(
+      cv = as.integer(cv),
+      random_state = as.integer(random_state),
+      max_iter = 1000L
+    )
+  } else if (is.null(model_propensity)) {
+    model_prop_py <- NULL
+  } else {
+    model_prop_py <- model_propensity
+  }
+
+  # Regression model
+  if (identical(model_regression, "auto")) {
+    RidgeCV <- sklearn$linear_model$RidgeCV
+    model_reg_py <- RidgeCV(cv = as.integer(cv))
+  } else if (is.null(model_regression)) {
+    model_reg_py <- NULL
+  } else {
+    model_reg_py <- model_regression
+  }
+
+  # Final model
+  if (identical(model_final, "auto")) {
+    LassoCV <- sklearn$linear_model$LassoCV
+    model_final_py <- LassoCV(
+      cv = as.integer(cv),
+      random_state = as.integer(random_state)
+    )
+  } else if (is.null(model_final)) {
+    model_final_py <- NULL
+  } else {
+    model_final_py <- model_final
+  }
+
   X_py <- reticulate::r_to_py(as.matrix(X))
   T_py <- reticulate::r_to_py(as.numeric(T))
   Y_py <- reticulate::r_to_py(as.numeric(Y))
```

---

### 変更4: DRLearner コンストラクタに引数を渡す（FIX #1, #2, #3, #5）

```diff
-  dr <- DRLearner(cv = as.integer(cv))
+  dr <- DRLearner(
+    model_propensity = model_prop_py,
+    model_regression = model_reg_py,
+    model_final = model_final_py,
+    cv = as.integer(cv),
+    min_propensity = min_propensity,
+    random_state = as.integer(random_state)
+  )
```

---

### 変更5: ITE の標準化を元に戻す（FIX #4）

```diff
   ite <- reticulate::py_to_r(ite_py)
+
+  # 標準化を元に戻す
+  if (standardize_y) {
+    ite <- ite * Y_sd
+  }

   ate <- mean(ite)

-  y0_hat <- Y - T * ite
-  y1_hat <- Y + (1 - T) * ite
+  y0_hat <- Y_orig - T * ite
+  y1_hat <- Y_orig + (1 - T) * ite
```

---

### 変更6: メタデータに標準化情報を追加

```diff
   list(
     model = dr,
     res = list(
       ite = ite,
       ate = ate,
       y0_hat = y0_hat,
       y1_hat = y1_hat,
       summary = summary_stats,
       samples = list(
         y0 = y0_hat,
         y1 = y1_hat,
         ite = ite
       ),
-      metadata = NULL
+      metadata = list(
+        standardize_y = standardize_y,
+        Y_mean = Y_mean,
+        Y_sd = Y_sd,
+        min_propensity = min_propensity
+      )
     )
   )
```

---

### 変更7: predict 関数で標準化を元に戻す

```diff
   dr <- object$fit$model
   ite_py <- dr$effect(X_new, T0 = 0L, T1 = 1L)
   ite <- reticulate::py_to_r(ite_py)

+  # 標準化を元に戻す
+  metadata <- object$fit$res$metadata
+  if (!is.null(metadata) && metadata$standardize_y) {
+    ite <- ite * metadata$Y_sd
+  }
+
   switch(
     type,
     "ate" = mean(ite),
     "ite" = ite,
     stop(sprintf("type '%s' not supported for newdata prediction", type))
   )
```

---

## 修正の影響範囲

### 破壊的変更
- **なし**（後方互換性を維持）
- デフォルト値を使えば既存コードはそのまま動作

### 新機能
- `min_propensity`: 傾向スコアクリッピング（デフォルト 0.05）
- `standardize_y`: Y の標準化（デフォルト TRUE）
- `model_final`: 最終段モデルの指定（デフォルト "auto" = LassoCV）

### パフォーマンス改善
- n=500, SexGene 条件での RMSE: **19.94 → <1.0** (予想)
- クリップ率: **100% → 0%** (予想)
- 全体的な安定性向上

---

## テスト方法

```bash
# 1. 修正版ファイルをパッケージに適用
cp sandbox/demov04/methods_drlearner_FIXED.R \
   packages/cfomics/R/methods_drlearner.R

# 2. パッケージ再ビルド
Rscript -e 'devtools::document("packages/cfomics")'
Rscript -e 'devtools::load_all("packages/cfomics")'

# 3. 検証スクリプト実行
Rscript sandbox/demov04/test_fix_n500.R

# 4. ベンチマーク全体を再実行
Rscript sandbox/demov04/run_cate_demov04.R
```

---

## 期待される結果

### 修正前（現状）
```
n=500, linear, SexGene, drlearner: RMSE = 19.94
  - 100% clipped (60 at -20, 40 at +20)
  - Complete failure
```

### 修正後（期待）
```
n=500, linear, SexGene, drlearner: RMSE < 1.0
  - 0% clipped
  - Normal prediction
  - 61-76% improvement vs SexOnly (matching other conditions)
```
