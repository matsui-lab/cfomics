# PR-2: Core Registry + as_cf_data Generic

## 目的

Method 拡張を registry ベースで行えるようにし、`as_cf_data()` を S3 generic として定義することで、将来のパッケージ分離（cfomicsPython, cfomicsAdapters 等）に備える。

## 前提条件

- [ ] PR-1 がマージ済み（packages/ 構造）
- [ ] main ブランチが安定

## ブランチ

```bash
git checkout main
git pull origin main
git checkout -b feature/core-registry
```

---

## 作業手順

### Step 1: registry.R の作成

**ファイル**: `packages/cfomics/R/registry.R`

```r
#' @title Method Registry for cfomics
#' @description Internal registry for causal inference methods
#' @name registry
#' @keywords internal

# Internal environment for registry storage
.cfomics_registry <- new.env(parent = emptyenv())
.cfomics_registry$methods <- list()

#' Register a causal inference method
#'
#' @param id Character. Method identifier (e.g., "grf", "dowhy_gcm")
#' @param fit_fun Function. The fitting function
#' @param predict_fun Function or NULL. The prediction function
#' @param requires_python Logical. Whether Python is required
#' @param package Character or NULL. Package providing this method
#' @param description Character or NULL. Human-readable description
#'
#' @return Invisible TRUE on success
#' @export
#'
#' @examples
#' \dontrun{
#' cf_register_method(
#'   id = "my_method",
#'   fit_fun = my_fit_function,
#'   predict_fun = my_predict_function,
#'   requires_python = FALSE,
#'   package = "mypackage",
#'   description = "My custom causal inference method"
#' )
#' }
cf_register_method <- function(id,
                                fit_fun,
                                predict_fun = NULL,
                                requires_python = FALSE,
                                package = NULL,
                                description = NULL) {
  # Validation

  if (!is.character(id) || length(id) != 1 || nchar(id) == 0) {
    rlang::abort(
      "Method `id` must be a non-empty character string",
      class = "cfomics_registry_error"
    )
  }

  if (!is.function(fit_fun)) {
    rlang::abort(
      "Method `fit_fun` must be a function",
      class = "cfomics_registry_error"
    )
  }

  if (!is.null(predict_fun) && !is.function(predict_fun)) {
    rlang::abort(
      "Method `predict_fun` must be a function or NULL",
      class = "cfomics_registry_error"
    )
  }

  # Register

.cfomics_registry$methods[[id]] <- list(
    id = id,
    fit_fun = fit_fun,
    predict_fun = predict_fun,
    requires_python = isTRUE(requires_python),
    package = package %||% "cfomics",
    description = description %||% "",
    registered_at = Sys.time()
  )

  invisible(TRUE)
}

#' List available causal inference methods
#'
#' @param include_unavailable Logical. Include methods with unavailable dependencies?
#'
#' @return A data.frame with method information
#' @export
#'
#' @examples
#' cf_methods()
#' cf_methods(include_unavailable = TRUE)
cf_methods <- function(include_unavailable = FALSE) {
  methods <- .cfomics_registry$methods

  if (length(methods) == 0) {
    return(data.frame(
      id = character(),
      package = character(),
      requires_python = logical(),
      available = logical(),
      description = character(),
      stringsAsFactors = FALSE
    ))
  }

  df <- do.call(rbind, lapply(methods, function(m) {
    # Check availability
    available <- if (m$requires_python) {
      # Check if Python and required modules are available
      tryCatch({
        cf_has_python(m$id)
      }, error = function(e) FALSE)
    } else {
      TRUE
    }

    data.frame(
      id = m$id,
      package = m$package,
      requires_python = m$requires_python,
      available = available,
      description = m$description,
      stringsAsFactors = FALSE
    )
  }))

  rownames(df) <- NULL

  if (!include_unavailable) {
    df <- df[df$available, , drop = FALSE]
  }

  df
}

#' Get a registered method by ID
#'
#' @param id Character. Method identifier
#'
#' @return Method list or NULL if not found
#' @keywords internal
#' @noRd
.cf_get_method <- function(id) {
  .cfomics_registry$methods[[id]]
}

#' Check if a method is registered
#'
#' @param id Character. Method identifier
#'
#' @return Logical
#' @keywords internal
#' @noRd
.cf_has_method <- function(id) {
  id %in% names(.cfomics_registry$methods)
}

#' Clear all registered methods (for testing)
#'
#' @keywords internal
#' @noRd
.cf_clear_registry <- function() {
  .cfomics_registry$methods <- list()
  invisible(TRUE)
}
```

### Step 2: as_cf_data.R の作成

**ファイル**: `packages/cfomics/R/as_cf_data.R`

```r
#' @title Convert to cf_data
#' @description Generic function to convert various data types to cf_data format
#'
#' @param x Input data (data.frame, SummarizedExperiment, etc.)
#' @param outcome_name Character. Name of the outcome variable
#' @param treatment_name Character. Name of the treatment variable
#' @param covariate_names Character vector. Names of covariate variables
#' @param ... Additional arguments passed to methods
#'
#' @return A cf_data object
#' @export
#'
#' @examples
#' df <- data.frame(
#'   Y = rnorm(100),
#'   T = rbinom(100, 1, 0.5),
#'   X1 = rnorm(100),
#'   X2 = rnorm(100)
#' )
#' cf_data <- as_cf_data(df, "Y", "T", c("X1", "X2"))
#' print(cf_data)
as_cf_data <- function(x, outcome_name, treatment_name, covariate_names, ...) {
  UseMethod("as_cf_data")
}

#' @rdname as_cf_data
#' @export
as_cf_data.data.frame <- function(x,
                                   outcome_name,
                                   treatment_name,
                                   covariate_names,
                                   ...) {
  # Validate inputs
  if (!is.character(outcome_name) || length(outcome_name) != 1) {
    rlang::abort(
      "`outcome_name` must be a single character string",
      class = "cfomics_validation_error"
    )
  }

  if (!is.character(treatment_name) || length(treatment_name) != 1) {
    rlang::abort(
      "`treatment_name` must be a single character string",
      class = "cfomics_validation_error"
    )
  }

  if (!is.character(covariate_names) || length(covariate_names) == 0) {
    rlang::abort(
      "`covariate_names` must be a non-empty character vector",
      class = "cfomics_validation_error"
    )
  }

  # Check columns exist
  all_vars <- c(outcome_name, treatment_name, covariate_names)
  missing_vars <- setdiff(all_vars, names(x))

  if (length(missing_vars) > 0) {
    rlang::abort(
      paste0("Variables not found in data: ", paste(missing_vars, collapse = ", ")),
      class = "cfomics_validation_error"
    )
  }

  # Check treatment is binary
  t_vals <- unique(x[[treatment_name]])
  t_vals <- t_vals[!is.na(t_vals)]
  if (!all(t_vals %in% c(0, 1))) {
    rlang::warn(
      paste0("Treatment variable '", treatment_name, "' is not binary (0/1). ",
             "Some methods may not work correctly."),
      class = "cfomics_treatment_warning"
    )
  }

  # Build cf_data object
  structure(
    list(
      data = x,
      outcome_name = outcome_name,
      treatment_name = treatment_name,
      covariate_names = covariate_names,
      meta = list(
        n = nrow(x),
        p = length(covariate_names),
        types = list(
          outcome = class(x[[outcome_name]])[1],
          treatment = class(x[[treatment_name]])[1],
          covariates = vapply(
            x[, covariate_names, drop = FALSE],
            function(col) class(col)[1],
            character(1)
          )
        ),
        n_missing = sum(!complete.cases(x[, all_vars])),
        treatment_prevalence = mean(x[[treatment_name]], na.rm = TRUE),
        created_at = Sys.time()
      )
    ),
    class = "cf_data"
  )
}

#' @rdname as_cf_data
#' @export
as_cf_data.default <- function(x, outcome_name, treatment_name, covariate_names, ...) {
  # Try to coerce to data.frame
  if (is.matrix(x)) {
    x <- as.data.frame(x)
    return(as_cf_data.data.frame(x, outcome_name, treatment_name, covariate_names, ...))
  }

  # Check if Bioconductor object (defer to cfomicsAdapters)
  if (inherits(x, "SummarizedExperiment") ||
      inherits(x, "MultiAssayExperiment") ||
      inherits(x, "SingleCellExperiment")) {
    rlang::abort(
      paste0(
        "Bioconductor objects require the 'cfomicsAdapters' package.\n",
        "Install with: remotes::install_github('matsui-lab/cfomics', subdir='packages/cfomicsAdapters')"
      ),
      class = "cfomics_missing_package_error"
    )
  }

  rlang::abort(
    paste0("Cannot convert object of class '", class(x)[1], "' to cf_data.\n",
           "Supported types: data.frame, matrix"),
    class = "cfomics_validation_error"
  )
}

#' @export
print.cf_data <- function(x, ...) {
  cli::cli_h1("cf_data object")
  cli::cli_ul()
  cli::cli_li("Samples: {.val {x$meta$n}}")
  cli::cli_li("Covariates: {.val {x$meta$p}}")
  cli::cli_li("Outcome: {.field {x$outcome_name}} ({x$meta$types$outcome})")
  cli::cli_li("Treatment: {.field {x$treatment_name}} ({x$meta$types$treatment})")
  cli::cli_li("Treatment prevalence: {.val {round(x$meta$treatment_prevalence * 100, 1)}%}")
  if (x$meta$n_missing > 0) {
    cli::cli_li("Missing observations: {.val {x$meta$n_missing}}")
  }
  cli::cli_end()
  invisible(x)
}

#' @export
summary.cf_data <- function(object, ...) {
  cat("cf_data Summary\n")
  cat(strrep("-", 40), "\n")
  cat("Dimensions:", object$meta$n, "samples x", object$meta$p, "covariates\n")
  cat("Outcome:", object$outcome_name, "\n")
  cat("Treatment:", object$treatment_name, "\n")
  cat("Covariates:", paste(object$covariate_names, collapse = ", "), "\n")
  cat("Treatment prevalence:", round(object$meta$treatment_prevalence * 100, 1), "%\n")
  cat("Missing observations:", object$meta$n_missing, "\n")
  invisible(object)
}

#' Validate cf_data object
#'
#' @param x Object to validate
#' @param strict Logical. If TRUE, perform stricter validation
#'
#' @return TRUE if valid, otherwise throws an error
#' @export
#'
#' @examples
#' df <- data.frame(Y = 1:10, T = rep(0:1, 5), X1 = rnorm(10))
#' cf_data <- as_cf_data(df, "Y", "T", "X1")
#' validate_cf_data(cf_data)
validate_cf_data <- function(x, strict = FALSE) {
  if (!inherits(x, "cf_data")) {
    rlang::abort(
      "Object must be of class 'cf_data'",
      class = "cfomics_validation_error"
    )
  }

  # Check required fields
  required_fields <- c("data", "outcome_name", "treatment_name", "covariate_names", "meta")
  missing_fields <- setdiff(required_fields, names(x))

  if (length(missing_fields) > 0) {
    rlang::abort(
      paste0("Missing required fields: ", paste(missing_fields, collapse = ", ")),
      class = "cfomics_validation_error"
    )
  }

  # Check data is data.frame
  if (!is.data.frame(x$data)) {
    rlang::abort(
      "`data` must be a data.frame",
      class = "cfomics_validation_error"
    )
  }

  # Check meta fields
  required_meta <- c("n", "p")
  missing_meta <- setdiff(required_meta, names(x$meta))

  if (length(missing_meta) > 0) {
    rlang::abort(
      paste0("Missing required meta fields: ", paste(missing_meta, collapse = ", ")),
      class = "cfomics_validation_error"
    )
  }

  # Strict validation
  if (strict) {
    # Check dimensions match
    if (x$meta$n != nrow(x$data)) {
      rlang::abort(
        "meta$n does not match nrow(data)",
        class = "cfomics_validation_error"
      )
    }

    if (x$meta$p != length(x$covariate_names)) {
      rlang::abort(
        "meta$p does not match length of covariate_names",
        class = "cfomics_validation_error"
      )
    }
  }

  invisible(TRUE)
}

#' Extract data.frame from cf_data
#'
#' @param x cf_data object
#' @param include Character. Which variables to include: "all", "model" (outcome + treatment + covariates), "covariates"
#'
#' @return A data.frame
#' @export
cf_data_frame <- function(x, include = c("all", "model", "covariates")) {
  validate_cf_data(x)
  include <- match.arg(include)

  switch(include,
    all = x$data,
    model = x$data[, c(x$outcome_name, x$treatment_name, x$covariate_names), drop = FALSE],
    covariates = x$data[, x$covariate_names, drop = FALSE]
  )
}
```

### Step 3: aaa-onload.R の更新

R-native methods を registry に登録する処理を追加：

**ファイル**: `packages/cfomics/R/aaa-onload.R` に追加

```r
# Add to .onLoad function:

.onLoad <- function(libname, pkgname) {
 # Existing code (if any)...

  # Register R-native methods
  .register_core_methods()

  invisible()
}

#' Register core R-native methods
#' @keywords internal
#' @noRd
.register_core_methods <- function() {
  # GRF (Generalized Random Forest)
  cf_register_method(
    id = "grf",
    fit_fun = cf_fit_grf,
    predict_fun = predict_cf_grf,
    requires_python = FALSE,
    package = "cfomics",
    description = "Generalized Random Forest for heterogeneous treatment effects"
  )

  # IPW (Inverse Probability Weighting)
  cf_register_method(
    id = "ipw",
    fit_fun = cf_fit_ipw,
    predict_fun = predict_cf_ipw,
    requires_python = FALSE,
    package = "cfomics",
    description = "Inverse Probability Weighting estimator"
  )

  # G-formula
  cf_register_method(
    id = "gformula",
    fit_fun = cf_fit_gformula,
    predict_fun = predict_cf_gformula,
    requires_python = FALSE,
    package = "cfomics",
    description = "G-computation formula estimator"
  )

  invisible()
}
```

### Step 4: cf_fit.R の修正

Registry-based dispatch に変更：

**ファイル**: `packages/cfomics/R/cf_fit.R` の `cf_fit()` 関数を修正

```r
#' @title Fit a causal inference model
#' @description Unified interface for causal inference methods
#'
#' @param formula A formula of the form Y ~ T | X1 + X2 + ...
#' @param data A data.frame or cf_data object
#' @param method Character. The method to use. See `cf_methods()` for available methods.
#' @param ... Additional arguments passed to the method-specific fitting function
#'
#' @return A cfomics_result object
#' @export
#'
#' @examples
#' \donttest{
#' # Generate example data
#' set.seed(42)
#' n <- 200
#' df <- data.frame(
#'   X1 = rnorm(n),
#'   X2 = rnorm(n),
#'   T = rbinom(n, 1, 0.5)
#' )
#' df$Y <- 1 + 2 * df$T + 0.5 * df$X1 + rnorm(n)
#'
#' # Fit with GRF
#' result <- cf_fit(Y ~ T | X1 + X2, data = df, method = "grf")
#' summary(result)
#' }
cf_fit <- function(formula, data, method = "grf", ...) {
  # Parse formula
  parsed <- parse_cf_formula(formula)

  # Convert to cf_data if needed
  if (!inherits(data, "cf_data")) {
    cf_data_obj <- as_cf_data(
      data,
      outcome_name = parsed$outcome,
      treatment_name = parsed$treatment,
      covariate_names = parsed$covariates
    )
  } else {
    cf_data_obj <- data
  }

  # Look up method in registry
  method_info <- .cf_get_method(method)

  if (is.null(method_info)) {
    # Provide helpful error message
    available <- cf_methods(include_unavailable = TRUE)

    if (nrow(available) == 0) {
      rlang::abort(
        paste0(
          "Method '", method, "' is not registered.\n",
          "No methods are currently available. ",
          "This may indicate a package loading issue."
        ),
        class = "cfomics_method_error"
      )
    }

    # Check if method exists but is unavailable
    if (method %in% available$id) {
      info <- available[available$id == method, ]
      if (info$requires_python) {
        rlang::abort(
          paste0(
            "Method '", method, "' requires Python but Python is not available.\n",
            "Setup Python with: cf_install_python_env()\n",
            "Or install the cfomicsPython package for full Python support."
          ),
          class = "cfomics_python_error"
        )
      }
    }

    # Method not found at all
    available_methods <- paste(available$id, collapse = ", ")
    rlang::abort(
      paste0(
        "Method '", method, "' is not registered.\n",
        "Available methods: ", available_methods, "\n",
        "For Python methods, install cfomicsPython package."
      ),
      class = "cfomics_method_error"
    )
  }

  # Check Python availability if needed
  if (method_info$requires_python) {
    cf_require_python(method)
  }

  # Call the fit function
  result <- method_info$fit_fun(
    formula = formula,
    data = data,
    ...
  )

  # Validate result
  validate_cfomics_result(result)

  result
}
```

### Step 5: テストファイルの作成

**ファイル**: `packages/cfomics/tests/testthat/test-registry.R`

```r
test_that("cf_register_method works", {
  # Clear registry for testing
  .cf_clear_registry()

  # Register a test method
  result <- cf_register_method(
    id = "test_method",
    fit_fun = function(...) list(),
    requires_python = FALSE,
    description = "Test method"
  )

  expect_true(result)
  expect_true(.cf_has_method("test_method"))

  # Clean up
  .cf_clear_registry()
  .register_core_methods()
})

test_that("cf_register_method validates inputs", {
  expect_error(
    cf_register_method(id = "", fit_fun = identity),
    class = "cfomics_registry_error"
  )

  expect_error(
    cf_register_method(id = "test", fit_fun = "not a function"),
    class = "cfomics_registry_error"
  )
})

test_that("cf_methods returns registered methods", {
  methods <- cf_methods()

  expect_s3_class(methods, "data.frame")
  expect_true("id" %in% names(methods))
  expect_true("grf" %in% methods$id)
  expect_true("ipw" %in% methods$id)
  expect_true("gformula" %in% methods$id)
})

test_that("cf_methods respects include_unavailable", {
  # All methods should be available (R-native)
  methods_all <- cf_methods(include_unavailable = TRUE)
  methods_avail <- cf_methods(include_unavailable = FALSE)

  # R-native methods should be in both
  expect_true("grf" %in% methods_avail$id)
})

test_that(".cf_get_method returns correct method", {
  method <- .cf_get_method("grf")

  expect_type(method, "list")
  expect_equal(method$id, "grf")
  expect_true(is.function(method$fit_fun))
  expect_false(method$requires_python)
})

test_that(".cf_get_method returns NULL for unknown method", {
  method <- .cf_get_method("nonexistent_method")
  expect_null(method)
})
```

**ファイル**: `packages/cfomics/tests/testthat/test-cf_data.R`

```r
test_that("as_cf_data.data.frame creates valid cf_data", {
  df <- data.frame(
    Y = rnorm(100),
    T = rbinom(100, 1, 0.5),
    X1 = rnorm(100),
    X2 = rnorm(100)
  )

  result <- as_cf_data(df, "Y", "T", c("X1", "X2"))

  expect_s3_class(result, "cf_data")
  expect_equal(result$outcome_name, "Y")
  expect_equal(result$treatment_name, "T")
  expect_equal(result$covariate_names, c("X1", "X2"))
  expect_equal(result$meta$n, 100)
  expect_equal(result$meta$p, 2)
})

test_that("as_cf_data validates inputs", {
  df <- data.frame(Y = 1:10, T = 0:1, X1 = rnorm(10))

  # Missing column
 expect_error(
    as_cf_data(df, "Y", "T", c("X1", "X2")),
    class = "cfomics_validation_error"
  )

  # Invalid outcome_name
  expect_error(
    as_cf_data(df, c("Y", "T"), "T", "X1"),
    class = "cfomics_validation_error"
  )

  # Empty covariates
  expect_error(
    as_cf_data(df, "Y", "T", character(0)),
    class = "cfomics_validation_error"
  )
})

test_that("as_cf_data warns for non-binary treatment", {
  df <- data.frame(Y = 1:10, T = 1:10, X1 = rnorm(10))

  expect_warning(
    as_cf_data(df, "Y", "T", "X1"),
    class = "cfomics_treatment_warning"
  )
})

test_that("validate_cf_data works", {
  df <- data.frame(Y = 1:10, T = rep(0:1, 5), X1 = rnorm(10))
  cf_data <- as_cf_data(df, "Y", "T", "X1")

  expect_true(validate_cf_data(cf_data))
  expect_true(validate_cf_data(cf_data, strict = TRUE))
})

test_that("validate_cf_data catches invalid objects", {
  expect_error(
    validate_cf_data(list(data = data.frame())),
    class = "cfomics_validation_error"
  )
})

test_that("cf_data_frame extracts data correctly", {
  df <- data.frame(Y = 1:10, T = rep(0:1, 5), X1 = rnorm(10), X2 = rnorm(10), Extra = 1:10)
  cf_data <- as_cf_data(df, "Y", "T", c("X1", "X2"))

  # All columns
  result_all <- cf_data_frame(cf_data, "all")
  expect_equal(ncol(result_all), 5)

  # Model columns only
  result_model <- cf_data_frame(cf_data, "model")
  expect_equal(names(result_model), c("Y", "T", "X1", "X2"))

  # Covariates only
  result_cov <- cf_data_frame(cf_data, "covariates")
  expect_equal(names(result_cov), c("X1", "X2"))
})

test_that("print.cf_data works", {
  df <- data.frame(Y = 1:10, T = rep(0:1, 5), X1 = rnorm(10))
  cf_data <- as_cf_data(df, "Y", "T", "X1")

  expect_output(print(cf_data), "cf_data object")
})

test_that("as_cf_data.default handles matrices", {
  mat <- matrix(rnorm(40), ncol = 4)
  colnames(mat) <- c("Y", "T", "X1", "X2")

  result <- as_cf_data(mat, "Y", "T", c("X1", "X2"))
  expect_s3_class(result, "cf_data")
})

test_that("as_cf_data.default errors for unsupported types", {
  expect_error(
    as_cf_data(list(a = 1), "Y", "T", "X1"),
    class = "cfomics_validation_error"
  )
})
```

### Step 6: NAMESPACE の更新

```bash
cd packages/cfomics
Rscript -e "devtools::document()"
```

追加される export:
- `cf_register_method`
- `cf_methods`
- `as_cf_data`
- `validate_cf_data`
- `cf_data_frame`

### Step 7: 既存テストの修正

既存のテストが registry ベースの `cf_fit` で動くことを確認し、必要に応じて修正。

---

## 検証手順

### Step A: ドキュメント生成

```bash
cd packages/cfomics
Rscript -e "devtools::document()"
```

### Step B: テスト実行

```bash
cd packages/cfomics
Rscript -e "devtools::test()"
```

期待: 全テストパス

### Step C: パッケージチェック

```bash
cd packages/cfomics
R CMD build .
R CMD check cfomics_*.tar.gz --as-cran
```

### Step D: 動作確認

```r
library(cfomics)

# Registry の確認
cf_methods()

# cf_data の作成
df <- data.frame(Y = rnorm(100), T = rbinom(100, 1, 0.5), X1 = rnorm(100))
cf_data <- as_cf_data(df, "Y", "T", "X1")
print(cf_data)

# cf_fit (registry 経由)
result <- cf_fit(Y ~ T | X1, data = df, method = "gformula")
summary(result)
```

---

## DoD チェックリスト

- [ ] `registry.R` が作成されている
- [ ] `as_cf_data.R` が作成されている
- [ ] `cf_register_method()` が動作する
- [ ] `cf_methods()` が R-native 3 種を返す
- [ ] `as_cf_data()` が data.frame を変換できる
- [ ] `cf_fit()` が registry 経由で dispatch する
- [ ] 既存のテストが全てパスする
- [ ] 新規テスト（test-registry.R, test-cf_data.R）がパスする
- [ ] `R CMD check` が ERROR なしで通る

---

## コミット手順

### Commit 1: registry.R

```bash
git add packages/cfomics/R/registry.R
git commit -m "feat(core): add method registry system

- Add cf_register_method() for dynamic method registration
- Add cf_methods() to list available methods
- Add internal helpers for registry management

🤖 Generated with [Claude Code](https://claude.com/claude-code)

Co-Authored-By: Claude <noreply@anthropic.com>"
```

### Commit 2: as_cf_data.R

```bash
git add packages/cfomics/R/as_cf_data.R
git commit -m "feat(core): add as_cf_data generic and cf_data class

- Add as_cf_data() generic for input standardization
- Add as_cf_data.data.frame() method
- Add validate_cf_data() for contract validation
- Add cf_data_frame() for data extraction

🤖 Generated with [Claude Code](https://claude.com/claude-code)

Co-Authored-By: Claude <noreply@anthropic.com>"
```

### Commit 3: cf_fit refactor

```bash
git add packages/cfomics/R/cf_fit.R packages/cfomics/R/aaa-onload.R
git commit -m "refactor(core): migrate cf_fit to registry-based dispatch

- Update cf_fit() to use method registry
- Register R-native methods (grf, ipw, gformula) in .onLoad()
- Improve error messages for missing methods

🤖 Generated with [Claude Code](https://claude.com/claude-code)

Co-Authored-By: Claude <noreply@anthropic.com>"
```

### Commit 4: Tests

```bash
git add packages/cfomics/tests/testthat/test-registry.R packages/cfomics/tests/testthat/test-cf_data.R
git commit -m "test(core): add tests for registry and cf_data

- Add test-registry.R for method registration
- Add test-cf_data.R for cf_data creation/validation

🤖 Generated with [Claude Code](https://claude.com/claude-code)

Co-Authored-By: Claude <noreply@anthropic.com>"
```

### Commit 5: Documentation

```bash
git add packages/cfomics/NAMESPACE packages/cfomics/man/
git commit -m "docs: regenerate documentation for new exports

🤖 Generated with [Claude Code](https://claude.com/claude-code)

Co-Authored-By: Claude <noreply@anthropic.com>"
```

---

## PR 作成

```bash
git push -u origin feature/core-registry

gh pr create --title "feat(core): add method registry and cf_data contract" --body "$(cat <<'EOF'
## Summary

- Add method registry system for dynamic method registration
- Add `cf_data` class as standardized input format
- Refactor `cf_fit()` to use registry-based dispatch

## New APIs

### Method Registry
- `cf_register_method()` - Register a new causal inference method
- `cf_methods()` - List available methods

### cf_data
- `as_cf_data()` - Convert data to standardized format
- `validate_cf_data()` - Validate cf_data objects
- `cf_data_frame()` - Extract data.frame from cf_data

## Breaking Changes

None. Existing `cf_fit()` calls work unchanged.

## Test plan

- [ ] `devtools::test()` passes
- [ ] `R CMD check` passes
- [ ] `cf_methods()` returns grf, ipw, gformula
- [ ] `cf_fit()` works with all R-native methods

🤖 Generated with [Claude Code](https://claude.com/claude-code)
EOF
)"
```

---

## 次のステップ

PR-2 がマージされたら:
- **PR-3**: cfomicsPython skeleton の作成
- **PR-4**: Python methods の移植

---

*Last updated: 2026-01-05*
