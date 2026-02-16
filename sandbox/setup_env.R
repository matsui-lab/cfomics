# setup_env.R
# cfomics ベンチマーク環境のセットアップスクリプト
# 実行: source("sandbox/setup_env.R")

# Python環境の指定（sandbox/.venv を使用）
VENV_PYTHON <- file.path(dirname(rstudioapi::getActiveProject() %||%
                           normalizePath(".")), "sandbox", ".venv", "bin", "python3")

# フォールバック: 相対パスで指定
if (!file.exists(VENV_PYTHON)) {
  VENV_PYTHON <- normalizePath("sandbox/.venv/bin/python3", mustWork = FALSE)
}

Sys.setenv(RETICULATE_PYTHON = VENV_PYTHON)
Sys.setenv(TF_CPP_MIN_LOG_LEVEL = "3")  # TensorFlowの冗長ログを抑制

# cfomics をロード
devtools::load_all("packages/cfomics", quiet = TRUE)

cat("=== cfomics ベンチマーク環境 ===\n")
cat("Python:", VENV_PYTHON, "\n")
cat("利用可能メソッド: GRF, DRLearner, GANITE, CAVAE(CEVAE)\n")
cat("モジュール確認:\n")
cat("  econml (DRLearner):", reticulate::py_module_available("econml"), "\n")
cat("  torch  (CEVAE)    :", reticulate::py_module_available("torch"), "\n")
cat("  pyro   (CEVAE)    :", reticulate::py_module_available("pyro"), "\n")
cat("  tensorflow (GANITE):", reticulate::py_module_available("tensorflow"), "\n")
