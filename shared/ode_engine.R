# ============================================================================
# mrgsolve model files, solved in base R.
#
# Every model in this repository is written once, as an mrgsolve model file
# (models/*.cpp). Run locally, mrgsolve compiles and solves it. A browser
# can't compile C++, so for the web version this file reads the same model
# file, translates its C++ into vectorised R, and integrates it with a stiff
# Rosenbrock solver. tests/ checks the two against each other.
#
# Supported subset of the mrgsolve model syntax:
#   blocks   $PARAM $CMT $INIT $MAIN $ODE $TABLE $CAPTURE $GLOBAL
#            (also written [PARAM] etc.), with or without @annotated
#   $GLOBAL  object-like #define macros
#   code     double x = expr;   x = expr;   dxdt_CMT = expr;   CMT_0 = expr;
#            F_CMT = expr;   if (cond) x = expr;   if (cond) { ... } else { ... }
#            cond ? a : b as the whole right-hand side of an assignment
#   maths    + - * / pow exp log log10 sqrt fabs fmax fmin, && || !
#   time     SOLVERTIME (or TIME) in $ODE
#
# Everything is vectorised over "columns": one column per subject (and,
# inside the solver, per finite-difference perturbation), so a population or
# a Jacobian is one call.
# ============================================================================

# --- Parsing ------------------------------------------------------------------

# Remove // and /* */ comments the way a C++ compiler does: whichever opens
# first wins, so "//* note */" is a line comment, not a block.
.mrg_strip_comments <- function(lines) {
  txt <- paste(lines, collapse = "\n")
  out <- character(0)
  i <- 1
  n <- nchar(txt)
  repeat {
    rest <- substr(txt, i, n)
    m <- regexpr("//|/\\*", rest)
    if (m < 0) { out <- c(out, rest); break }
    out <- c(out, substr(rest, 1, m - 1))
    if (substr(rest, m, m + 1) == "//") {
      nl <- regexpr("\n", substr(rest, m, nchar(rest)), fixed = TRUE)
      if (nl < 0) break
      i <- i + m - 1 + nl - 1           # keep the newline
    } else {
      close <- regexpr("*/", substr(rest, m + 2, nchar(rest)), fixed = TRUE)
      if (close < 0) break
      skipped <- substr(rest, m, m + 2 + close)
      out <- c(out, strrep("\n", lengths(regmatches(skipped, gregexpr("\n", skipped)))))
      i <- i + m + 2 + close + 1 - 1
    }
  }
  strsplit(paste(out, collapse = ""), "\n", fixed = TRUE)[[1]]
}

.mrg_blocks <- function(lines) {
  hdr <- grepl("^\\s*(\\$[A-Za-z]+|\\[\\s*[A-Za-z]+\\s*\\])", lines)
  blocks <- list()
  cur <- NULL
  for (i in seq_along(lines)) {
    if (hdr[i]) {
      name <- toupper(sub("^\\s*\\$?\\[?\\s*([A-Za-z]+).*$", "\\1", lines[i]))
      rest <- sub("^\\s*(\\$[A-Za-z]+|\\[\\s*[A-Za-z]+\\s*\\])", "", lines[i])
      cur <- list(name = name, annotated = grepl("@annotated", rest),
                  lines = sub("@[A-Za-z]+", "", rest))
      blocks[[length(blocks) + 1]] <- cur
    } else if (length(blocks)) {
      blocks[[length(blocks)]]$lines <- c(blocks[[length(blocks)]]$lines, lines[i])
    }
  }
  blocks
}

.mrg_eval_num <- function(expr, env) {
  as.numeric(eval(parse(text = .mrg_translate_expr(expr)), envir = env))
}

# NAME = value pairs (comma or newline separated) or annotated NAME : value : note
.mrg_pairs <- function(block, env = new.env()) {
  out <- list()
  body <- block$lines
  if (block$annotated) {
    for (l in body) {
      if (!nzchar(trimws(l)) || !grepl(":", l, fixed = TRUE)) next
      parts <- strsplit(l, ":", fixed = TRUE)[[1]]
      nm <- trimws(parts[1])
      out[[nm]] <- .mrg_eval_num(trimws(parts[2]), env)
      assign(nm, out[[nm]], envir = env)
    }
  } else {
    items <- unlist(strsplit(paste(body, collapse = ","), "[,\n]"))
    for (it in items) {
      it <- trimws(it)
      if (!nzchar(it)) next
      kv <- strsplit(it, "=", fixed = TRUE)[[1]]
      nm <- trimws(kv[1])
      out[[nm]] <- .mrg_eval_num(trimws(paste(kv[-1], collapse = "=")), env)
      assign(nm, out[[nm]], envir = env)
    }
  }
  unlist(out)
}

.mrg_names <- function(block) {
  if (block$annotated) {
    nm <- vapply(block$lines, function(l) trimws(strsplit(l, ":", fixed = TRUE)[[1]][1]), "")
    nm <- unname(nm[!is.na(nm) & nzchar(nm)])
  } else {
    nm <- unlist(strsplit(paste(block$lines, collapse = " "), "[[:space:],]+"))
    nm <- nm[nzchar(nm)]
  }
  nm
}

# Object-like macros (#define NAME body) are returned as strings; function-
# like ones (#define NAME(a, b) body) as list(args, body) under attr "fn".
# max/min macros are left alone: they are translated to pmax/pmin.
.mrg_defines <- function(block) {
  defs <- list()
  fn <- list()
  for (l in block$lines) {
    m <- regmatches(l, regexec("^\\s*#define\\s+([A-Za-z_][A-Za-z0-9_]*)\\s+(.+)$", l))[[1]]
    if (length(m) == 3) { defs[[m[2]]] <- trimws(m[3]); next }
    f <- regmatches(l, regexec("^\\s*#define\\s+([A-Za-z_][A-Za-z0-9_]*)\\(([^)]*)\\)\\s+(.+)$", l))[[1]]
    if (length(f) == 4 && !f[2] %in% c("max", "min")) {
      fn[[f[2]]] <- list(args = trimws(strsplit(f[3], ",")[[1]]), body = trimws(f[4]))
    }
  }
  if (length(fn)) defs[[".fn"]] <- fn
  defs
}

# Expand function-like macros at their call sites, as the preprocessor does.
.mrg_expand_fn <- function(code, fn) {
  if (!length(fn)) return(code)
  vapply(code, function(line) {
    for (k in 1:10) {
      hit <- FALSE
      for (nm in names(fn)) {
        pos <- regexpr(sprintf("\\b%s\\s*\\(", nm), line, perl = TRUE)
        if (pos < 0) next
        hit <- TRUE
        open <- pos + attr(pos, "match.length") - 1
        chars <- strsplit(line, "", fixed = TRUE)[[1]]
        depth <- 0; args <- character(0); start <- open + 1; j <- open
        repeat {
          c <- chars[j]
          if (c == "(") depth <- depth + 1
          if (c == ")") { depth <- depth - 1; if (depth == 0) break }
          if (c == "," && depth == 1) { args <- c(args, substr(line, start, j - 1)); start <- j + 1 }
          j <- j + 1
        }
        args <- c(args, substr(line, start, j - 1))
        body <- fn[[nm]]$body
        for (a in seq_along(fn[[nm]]$args)) {
          body <- gsub(sprintf("\\b%s\\b", fn[[nm]]$args[a]), paste0("(", trimws(args[a]), ")"), body, perl = TRUE)
        }
        line <- paste0(substr(line, 1, pos - 1), "(", body, ")", substr(line, j + 1, nchar(line)))
      }
      if (!hit) break
    }
    line
  }, "", USE.NAMES = FALSE)
}

.mrg_apply_defines <- function(code, defs) {
  if (!length(defs)) return(code)
  code <- .mrg_expand_fn(code, defs[[".fn"]])
  defs[[".fn"]] <- NULL
  for (k in 1:5) {                      # macros may reference macros
    before <- code
    for (nm in names(defs)) {
      code <- gsub(sprintf("\\b%s\\b", nm), paste0("(", defs[[nm]], ")"), code, perl = TRUE)
    }
    if (identical(code, before)) break
  }
  code
}

# --- Translation of C expressions and statements --------------------------------

# pow(a, b) -> ((a)^(b)): R's ^ is a primitive, a pow() closure is not, and
# the difference dominates the cost of models built from Hill functions.
.mrg_replace_pow <- function(x) {
  repeat {
    pos <- regexpr("\\bpow\\s*\\(", x, perl = TRUE)
    if (pos < 0) return(x)
    open <- pos + attr(pos, "match.length") - 1
    chars <- strsplit(x, "", fixed = TRUE)[[1]]
    depth <- 0; comma <- NA; j <- open
    repeat {
      c <- chars[j]
      if (c == "(") depth <- depth + 1
      if (c == ")") { depth <- depth - 1; if (depth == 0) break }
      if (c == "," && depth == 1 && is.na(comma)) comma <- j
      j <- j + 1
    }
    a <- substr(x, open + 1, comma - 1)
    b <- substr(x, comma + 1, j - 1)
    x <- paste0(substr(x, 1, pos - 1), "((", a, ")^(", b, "))", substr(x, j + 1, nchar(x)))
  }
}

.mrg_translate_expr <- function(x) {
  x <- gsub("std::", "", x, fixed = TRUE)
  x <- .mrg_replace_pow(x)
  x <- gsub("&&", "&", x, fixed = TRUE)
  x <- gsub("||", "|", x, fixed = TRUE)
  x <- gsub("\\bfabs\\(", "abs(", x, perl = TRUE)
  # pmax.int / pmin.int: internal, much cheaper than the pmax closure
  x <- gsub("\\bfmax\\(", "pmax.int(", x, perl = TRUE)
  x <- gsub("\\bfmin\\(", "pmin.int(", x, perl = TRUE)
  x <- gsub("\\bmax\\(", "pmax.int(", x, perl = TRUE)
  x <- gsub("\\bmin\\(", "pmin.int(", x, perl = TRUE)
  x <- gsub("\\bSOLVERTIME\\b", "t", x, perl = TRUE)
  x <- gsub("\\bTIME\\b", "t", x, perl = TRUE)
  x <- gsub("(\\d)[fF]\\b", "\\1", x, perl = TRUE)     # 1.0f literals
  # cond ? a : b, at the top level of the expression
  q <- .mrg_top_level(x, "?")
  if (!is.na(q)) {
    rest <- substring(x, q + 1)
    cpos <- .mrg_top_level(rest, ":")
    cond <- substring(x, 1, q - 1)
    a <- substring(rest, 1, cpos - 1)
    b <- substring(rest, cpos + 1)
    return(sprintf("ifelse(%s, %s, %s)", .mrg_translate_expr(cond),
                   .mrg_translate_expr(a), .mrg_translate_expr(b)))
  }
  trimws(x)
}

# Position of the first character ch at parenthesis depth 0, or NA.
.mrg_top_level <- function(x, ch) {
  chars <- strsplit(x, "", fixed = TRUE)[[1]]
  depth <- 0
  for (i in seq_along(chars)) {
    c <- chars[i]
    if (c == "(") depth <- depth + 1
    else if (c == ")") depth <- depth - 1
    else if (c == ch && depth == 0) return(i)
  }
  NA_integer_
}

# Split C code into statements, keeping if / else structure.
.mrg_statements <- function(code) {
  s <- paste(code, collapse = " ")
  out <- list()
  i <- 1
  n <- nchar(s)
  read_until_semicolon <- function(i) {
    j <- i
    depth <- 0
    while (j <= n) {
      c <- substr(s, j, j)
      if (c == "(") depth <- depth + 1
      if (c == ")") depth <- depth - 1
      if (c == ";" && depth == 0) break
      j <- j + 1
    }
    list(stmt = trimws(substr(s, i, j - 1)), next_i = j + 1)
  }
  read_parens <- function(i) {           # s[i] == "("
    depth <- 0
    j <- i
    repeat {
      c <- substr(s, j, j)
      if (c == "(") depth <- depth + 1
      if (c == ")") { depth <- depth - 1; if (depth == 0) break }
      j <- j + 1
    }
    list(inner = substr(s, i + 1, j - 1), next_i = j + 1)
  }
  read_body <- function(i) {             # a { block } or a single statement
    while (substr(s, i, i) %in% c(" ", "\t")) i <- i + 1
    if (substr(s, i, i) == "{") {
      depth <- 0
      j <- i
      repeat {
        c <- substr(s, j, j)
        if (c == "{") depth <- depth + 1
        if (c == "}") { depth <- depth - 1; if (depth == 0) break }
        j <- j + 1
      }
      list(stmts = .mrg_statements(substr(s, i + 1, j - 1)), next_i = j + 1)
    } else {
      r <- read_until_semicolon(i)
      list(stmts = list(list(type = "stmt", text = r$stmt)), next_i = r$next_i)
    }
  }
  while (i <= n) {
    while (i <= n && substr(s, i, i) %in% c(" ", "\t", "\n", ";")) i <- i + 1
    if (i > n) break
    if (grepl("^if\\s*\\(", substr(s, i, n))) {
      p <- regexpr("\\(", substr(s, i, n)) + i - 1
      cond <- read_parens(p)
      body <- read_body(cond$next_i)
      i <- body$next_i
      else_body <- NULL
      rest <- substr(s, i, n)
      if (grepl("^\\s*else\\b", rest)) {
        i <- i + regexpr("else", rest) + 4
        else_body <- read_body(i)
        i <- else_body$next_i
        else_body <- else_body$stmts
      }
      out[[length(out) + 1]] <- list(type = "if", cond = cond$inner,
                                     then = body$stmts, otherwise = else_body)
    } else {
      r <- read_until_semicolon(i)
      if (nzchar(r$stmt)) out[[length(out) + 1]] <- list(type = "stmt", text = r$stmt)
      i <- r$next_i
    }
  }
  out
}

# Translate statements to R. Conditional assignments become ifelse() so the
# code stays vectorised over columns.
.mrg_to_r <- function(stmts, cond = NULL) {
  out <- character(0)
  for (st in stmts) {
    if (st$type == "if") {
      cnd <- .mrg_translate_expr(st$cond)
      c_then <- if (is.null(cond)) sprintf("(%s)", cnd) else sprintf("(%s) & (%s)", cond, cnd)
      out <- c(out, .mrg_to_r(st$then, c_then))
      if (!is.null(st$otherwise)) {
        c_else <- if (is.null(cond)) sprintf("!(%s)", cnd) else sprintf("(%s) & !(%s)", cond, cnd)
        out <- c(out, .mrg_to_r(st$otherwise, c_else))
      }
      next
    }
    txt <- sub("^(capture|double|float|int|bool|const double)\\s+", "", st$text)
    eq <- regexpr("(?<![=!<>])=(?!=)", txt, perl = TRUE)
    if (eq < 0) next
    lhs <- trimws(substr(txt, 1, eq - 1))
    rhs <- .mrg_translate_expr(substr(txt, eq + 1, nchar(txt)))
    if (is.null(cond)) {
      out <- c(out, sprintf("%s <- %s", lhs, rhs))
    } else {
      out <- c(out, sprintf("%s <- ifelse(%s, %s, if (exists(\"%s\", inherits = FALSE)) %s else 0)",
                            lhs, cond, rhs, lhs, lhs))
    }
  }
  out
}

pow <- function(a, b) a^b

#' Read an mrgsolve model file
#'
#' @return an object with the parameters, compartments, initial values and
#'   the translated R code of $MAIN, $ODE and $TABLE
mrg_read <- function(path) {
  lines <- .mrg_strip_comments(readLines(path, warn = FALSE))
  blocks <- .mrg_blocks(lines)
  get <- function(nm) Filter(function(b) b$name == nm, blocks)
  defs <- do.call(c, lapply(get("GLOBAL"), .mrg_defines))

  penv <- new.env()
  param <- do.call(c, lapply(get("PARAM"), .mrg_pairs, env = penv))
  cmt <- do.call(c, lapply(get("CMT"), .mrg_names))
  init <- do.call(c, lapply(get("INIT"), .mrg_pairs, env = penv))
  cmt <- c(cmt, setdiff(names(init), cmt))
  init_all <- stats::setNames(rep(0, length(cmt)), cmt)
  if (length(init)) init_all[names(init)] <- init

  declared <- character(0)
  code <- function(nm) {
    b <- get(nm)
    if (!length(b)) return(character(0))
    src <- unlist(lapply(b, `[[`, "lines"))
    # `capture x = ...;` both assigns and requests x as an output
    caps <- regmatches(src, gregexpr("\\bcapture\\s+[A-Za-z_][A-Za-z0-9_]*", src))
    declared <<- c(declared, sub("^capture\\s+", "", unlist(caps)))
    .mrg_to_r(.mrg_statements(.mrg_apply_defines(src, defs)))
  }
  capture <- unique(unlist(lapply(get("CAPTURE"), .mrg_names)))
  main <- code("MAIN"); ode <- code("ODE"); table <- code("TABLE")

  # A captured macro (a #define for a concentration, say) becomes table code.
  macro_caps <- intersect(capture, names(defs))
  if (length(macro_caps)) {
    table <- c(table, sprintf("%s <- %s", macro_caps, vapply(macro_caps, function(k) {
      .mrg_translate_expr(.mrg_apply_defines(defs[[k]], defs))
    }, "")))
  }

  structure(list(
    path = path, param = param, cmt = cmt, init = init_all,
    main = main, ode = ode, table = table,
    capture = unique(c(declared, capture))
  ), class = "mrg_r_model")
}

print.mrg_r_model <- function(x, ...) {
  cat(sprintf("mrgsolve model %s: %d parameters, %d compartments\n",
              basename(x$path), length(x$param), length(x$cmt)))
  invisible(x)
}

# --- Evaluation --------------------------------------------------------------------

# Run translated code in an environment holding parameters, states and
# earlier results; return that environment.
.mrg_run <- function(code, vars) {
  env <- list2env(vars, parent = globalenv())
  env$pow <- pow
  eval(parse(text = code), envir = env)
  env
}

#' Per-subject quantities from $MAIN: derived variables, initial values
#' and bioavailability
#'
#' @param P data frame of parameters, one row per subject (any parameter not
#'   given takes the model default)
mrg_main <- function(model, P) {
  m <- nrow(P)
  pars <- lapply(names(model$param), function(nm) {
    if (nm %in% names(P)) P[[nm]] else rep(model$param[[nm]], m)
  })
  names(pars) <- names(model$param)
  # Compartment initial values are visible to the code as CMT_0.
  for (k in seq_along(model$cmt)) pars[[paste0(model$cmt[k], "_0")]] <- rep(model$init[[k]], m)
  env <- .mrg_run(model$main, pars)
  vars <- mget(ls(env), envir = env)
  vars <- vars[vapply(vars, is.numeric, TRUE)]
  vars <- lapply(vars, rep_len, m)

  init <- matrix(rep(model$init, m), nrow = length(model$cmt),
                 dimnames = list(model$cmt, NULL))
  bio <- matrix(1, nrow = length(model$cmt), ncol = m, dimnames = list(model$cmt, NULL))
  for (i in seq_along(model$cmt)) {
    nm0 <- paste0(model$cmt[i], "_0")
    if (!is.null(vars[[nm0]])) init[i, ] <- vars[[nm0]]
    nmf <- paste0("F_", model$cmt[i])
    if (!is.null(vars[[nmf]])) bio[i, ] <- vars[[nmf]]
    vars[[nm0]] <- init[i, ]
  }
  list(vars = vars, init = init, bio = bio)
}

# Build a vectorised right-hand side: f(t, Y, V) with Y [cmt x columns] and V
# a list of per-column parameters and $MAIN variables. The generated function
# unpacks only the names its code uses and is byte-compiled, which matters:
# it is called a few times per solver step.
mrg_rhs <- function(model, vnames = names(model$param)) {
  cm <- model$cmt
  words <- unique(unlist(regmatches(model$ode, gregexpr("[A-Za-z_][A-Za-z0-9_.]*", model$ode))))
  used <- setdiff(intersect(words, vnames), cm)
  lines <- c(
    sprintf("%s <- V[[\"%s\"]]", used, used),
    sprintf("%s <- Y[%d, ]", cm, seq_along(cm)),
    model$ode,
    "out <- matrix(0, nrow(Y), ncol(Y))",
    unlist(lapply(seq_along(cm), function(i) {
      nm <- paste0("dxdt_", cm[i])
      if (any(grepl(sprintf("^%s <-", nm), model$ode))) sprintf("out[%d, ] <- %s", i, nm) else NULL
    })),
    "out"
  )
  f <- function(t, Y, V) NULL
  body(f) <- parse(text = paste(c("{", lines, "}"), collapse = "
"))[[1]]
  environment(f) <- list2env(list(pow = pow), parent = globalenv())
  compiler::cmpfun(f)
}

#' Does the right-hand side depend on time explicitly (SOLVERTIME)?
mrg_uses_time <- function(model) any(grepl("\\bt\\b", model$ode))

# --- Stiff solver ----------------------------------------------------------------------

#' One Rosenbrock run from t0 to t1 for every column at once
#'
#' The modified Rosenbrock (2,3) pair of Shampine & Reichelt (SIAM J Sci
#' Comput 1997;18:1), the method behind MATLAB's ode23s. L-stable, so it
#' handles the stiffness of binding reactions and fast turnover. The
#' Jacobian comes from finite differences, all columns and perturbations in
#' one call of the right-hand side. It is the expensive part, so it is kept
#' for up to `jac_every` accepted steps and refreshed at once after a
#' rejected step; the test suite checks the result against mrgsolve.
#'
#' @param uses_t whether the model's right-hand side depends on time
#'   explicitly (if not, the time-derivative term is zero and is skipped)
#' @param jac the Jacobian state returned by the previous segment, carried
#'   over when the state is continuous across the boundary (an output time,
#'   not a bolus), so that it is not recomputed at every output time
rosenbrock_segment <- function(f, Y, t0, t1, V, u, h, rtol, atol, max_steps = 2e5,
                               uses_t = TRUE, jac_every = 5, jac = NULL) {
  n <- nrow(Y); m <- ncol(Y)
  d <- 1 / (2 + sqrt(2)); e32 <- 6 + sqrt(2)
  t <- t0
  # Parameters repeated for the n + 1 finite-difference columns per subject.
  Vj <- lapply(V, function(v) if (length(v) == m) rep(v, each = n + 1) else v)
  fu <- function(t, Y, Vx, uu) f(t, Y, Vx) + uu
  uj <- u[, rep(seq_len(m), each = n + 1), drop = FALSE]
  steps <- 0
  I <- diag(n)
  Jl <- NULL
  jac_age <- Inf
  if (!is.null(jac)) { Jl <- jac$J; jac_age <- jac$age }
  while (t < t1 - 1e-12 * max(1, abs(t1))) {
    h <- min(h, t1 - t)
    if (jac_age >= jac_every) {
      del <- sqrt(.Machine$double.eps) * pmax(abs(Y), atol)
      Yp <- Y[, rep(seq_len(m), each = n + 1), drop = FALSE]
      for (s in seq_len(m)) {
        base <- (s - 1) * (n + 1)
        for (j in seq_len(n)) Yp[j, base + 1 + j] <- Yp[j, base + 1 + j] + del[j, s]
      }
      Fp <- fu(t, Yp, Vj, uj)
      F0 <- Fp[, seq(1, m * (n + 1), by = n + 1), drop = FALSE]
      Jl <- lapply(seq_len(m), function(s) {
        base <- (s - 1) * (n + 1)
        (Fp[, base + 1 + seq_len(n), drop = FALSE] - F0[, s]) / rep(del[, s], each = n)
      })
      jac_age <- 0
    } else {
      F0 <- fu(t, Y, V, u)
    }
    Tt <- if (uses_t) {
      dt <- sqrt(.Machine$double.eps) * max(abs(t), 1)
      (fu(t + dt, Y, V, u) - F0) / dt
    } else 0
    Winv <- lapply(Jl, function(J) tryCatch(solve(I - h * d * J, tol = 0), error = function(e) NULL))
    if (any(vapply(Winv, is.null, TRUE))) { h <- h / 4; jac_age <- Inf; next }
    mult <- function(Rhs) {
      out <- Rhs
      for (s in seq_len(m)) out[, s] <- Winv[[s]] %*% Rhs[, s]
      out
    }
    k1 <- mult(F0 + h * d * Tt)
    F1 <- fu(t + 0.5 * h, Y + 0.5 * h * k1, V, u)
    k2 <- mult(F1 - k1) + k1
    Ynew <- Y + h * k2
    F2 <- fu(t + h, Ynew, V, u)
    k3 <- mult(F2 - e32 * (k2 - F1) - 2 * (k1 - F0) + h * d * Tt)
    err <- h / 6 * (k1 - 2 * k2 + k3)
    sc <- atol + rtol * pmax(abs(Y), abs(Ynew))
    enorm <- max(sqrt(colMeans((err / sc)^2)))
    if (is.finite(enorm) && enorm <= 1) {
      t <- t + h
      Y <- Ynew
      jac_age <- jac_age + 1
    } else {
      jac_age <- Inf
    }
    fac <- if (is.finite(enorm)) 0.8 * max(enorm, 1e-10)^(-1 / 3) else 0.1
    h <- h * min(5, max(0.1, fac))
    steps <- steps + 1
    if (steps > max_steps) stop("ODE solver exceeded the maximum number of steps")
  }
  list(Y = Y, h = h, jac = list(J = Jl, age = jac_age))
}

#' Simulate a model
#'
#' @param model from mrg_read()
#' @param P data frame of parameters, one row per subject (may be a single
#'   row, or have zero columns to use the defaults)
#' @param events data frame of doses: time, cmt (name), amt, and optionally
#'   rate (0 = bolus), ii and addl for repeats, ID (row of P) to dose one
#'   subject only, and evid (1 = dose, the default; 8 = set the compartment
#'   to amt, as mrgsolve's replace event)
#' @param times output times
#' @param init optional starting state (named vector, or [cmt x subject] matrix)
#' @param nonneg evaluate the right-hand side at max(state, 0)
#' @return list(time, states [time x cmt x subject], table outputs as
#'   matrices [time x subject])
mrg_solve <- function(model, P = data.frame(row.names = 1), events = NULL, times,
                      rtol = 1e-6, atol = 1e-10, h0 = 1e-4, init = NULL,
                      nonneg = FALSE) {
  if (!nrow(P)) P <- data.frame(row.names = 1)
  m <- nrow(P)
  n <- length(model$cmt)
  mn <- mrg_main(model, P)
  V <- c(lapply(stats::setNames(names(model$param), names(model$param)), function(nm) {
    if (nm %in% names(P)) P[[nm]] else rep(model$param[[nm]], m)
  }), mn$vars)
  V <- V[!duplicated(names(V), fromLast = TRUE)]
  f <- mrg_rhs(model, names(V))
  if (nonneg) {
    # For networks of concentrations: evaluate the rates at max(y, 0), so a
    # tiny negative excursion inside a trial step can't flip the sign of a
    # saturable rate. The true solution never goes negative, so this only
    # changes rejected steps, not the answer.
    f0 <- f
    f <- function(t, Y, V) f0(t, pmax(Y, 0), V)
  }

  ev <- mrg_expand_events(events, m)
  times <- sort(unique(times))
  inf_on <- ev[ev$rate > 0, , drop = FALSE]
  grid <- sort(unique(c(times, ev$time, inf_on$time + inf_on$amt / inf_on$rate)))
  grid <- grid[grid >= min(times) & grid <= max(times)]
  if (min(grid) > 0 && nrow(ev) && any(ev$time < min(grid))) {
    grid <- sort(unique(c(0, ev$time[ev$time < min(grid)], grid)))
  }
  out_idx <- match(times, grid)

  Y <- mn$init
  # Optional starting state: a named vector (every subject) or a
  # [cmt x subject] matrix, overriding the model's initial values.
  if (!is.null(init)) {
    if (is.matrix(init)) {
      Y[rownames(init), ] <- init
    } else {
      Y[names(init), ] <- init
    }
  }
  states <- array(NA_real_, c(length(times), n, m), dimnames = list(NULL, model$cmt, NULL))
  h <- h0
  jac <- NULL
  for (j in seq_along(grid)) {
    tj <- grid[j]
    # bolus doses at tj
    b <- ev[abs(ev$time - tj) < 1e-10 & ev$rate == 0, , drop = FALSE]
    if (nrow(b)) jac <- NULL                  # the state jumps: a fresh Jacobian
    for (k in seq_len(nrow(b))) {
      ci <- match(b$cmt[k], model$cmt)
      cols <- if (is.na(b$ID[k])) seq_len(m) else b$ID[k]
      # evid 8 replaces the amount in the compartment, as in mrgsolve
      Y[ci, cols] <- if (b$evid[k] == 8) b$amt[k] else Y[ci, cols] + b$amt[k] * mn$bio[ci, cols]
    }
    o <- which(out_idx == j)
    if (length(o)) states[o, , ] <- Y
    if (j == length(grid)) break
    # infusion input over [tj, t_{j+1})
    u <- matrix(0, n, m)
    act <- inf_on[inf_on$time <= tj + 1e-10 & tj < inf_on$time + inf_on$amt / inf_on$rate - 1e-10, , drop = FALSE]
    for (k in seq_len(nrow(act))) {
      ci <- match(act$cmt[k], model$cmt)
      cols <- if (is.na(act$ID[k])) seq_len(m) else act$ID[k]
      u[ci, cols] <- u[ci, cols] + act$rate[k] * mn$bio[ci, cols]
    }
    # After a bolus into a state that drives fast dynamics the step size
    # adapts down on its own; restarting it from h0 at every dose only costs
    # steps, so keep it unless it is larger than the coming segment.
    seg <- rosenbrock_segment(f, Y, tj, grid[j + 1], V, u, h, rtol, atol,
                              uses_t = mrg_uses_time(model), jac = jac)
    Y <- seg$Y
    h <- seg$h
    jac <- seg$jac
  }
  res <- list(time = times, states = states)
  if (length(model$table) || length(model$capture)) {
    res <- c(res, mrg_table(model, states, times, V))
  }
  res
}

#' Doses expanded over ii / addl
mrg_expand_events <- function(events, m) {
  if (is.null(events) || !nrow(events)) {
    return(data.frame(time = numeric(0), cmt = character(0), amt = numeric(0),
                      rate = numeric(0), ID = integer(0), evid = numeric(0)))
  }
  e <- events
  if (is.null(e$rate)) e$rate <- 0
  if (is.null(e$evid)) e$evid <- 1
  if (is.null(e$ii)) e$ii <- 0
  if (is.null(e$addl)) e$addl <- 0
  if (is.null(e$ID)) e$ID <- NA_integer_
  rows <- lapply(seq_len(nrow(e)), function(i) {
    k <- 0:e$addl[i]
    data.frame(time = e$time[i] + k * e$ii[i], cmt = as.character(e$cmt[i]),
               amt = e$amt[i], rate = e$rate[i], ID = e$ID[i], evid = e$evid[i])
  })
  do.call(rbind, rows)
}

#' $TABLE and $CAPTURE outputs at every output time and subject
mrg_table <- function(model, states, times, V) {
  nt <- length(times); m <- dim(states)[3]
  vars <- lapply(V, function(v) rep(if (length(v) == m) v else rep_len(v, m), each = nt))
  for (i in seq_along(model$cmt)) vars[[model$cmt[i]]] <- as.vector(states[, i, ])
  vars$t <- rep(times, m)
  # $TABLE can use $ODE locals (as in mrgsolve), so evaluate $ODE first on
  # the output grid, then $TABLE.
  env <- .mrg_run(c(model$ode, model$table), vars)
  out <- list()
  for (nm in model$capture) {
    out[[nm]] <- matrix(rep_len(get(nm, envir = env), nt * m), nrow = nt)
  }
  out
}

# --- Stepping API -------------------------------------------------------------------

#' Set a model up for manual stepping
#'
#' For computations that decide as they go when to stop, such as a clotting
#' test that ends when the fibrin integral crosses a threshold.
mrg_prepare <- function(model, P = data.frame(row.names = 1), nonneg = FALSE) {
  if (!nrow(P)) P <- data.frame(row.names = 1)
  m <- nrow(P)
  mn <- mrg_main(model, P)
  V <- c(lapply(stats::setNames(names(model$param), names(model$param)), function(nm) {
    if (nm %in% names(P)) P[[nm]] else rep(model$param[[nm]], m)
  }), mn$vars)
  V <- V[!duplicated(names(V), fromLast = TRUE)]
  f <- mrg_rhs(model, names(V))
  if (nonneg) {
    f0 <- f
    f <- function(t, Y, V) f0(t, pmax(Y, 0), V)
  }
  Y <- mn$init
  rownames(Y) <- model$cmt
  list(model = model, f = f, V = V, init = Y, m = m)
}

#' Advance a prepared model from t0 to t1 (no dosing inside the interval)
mrg_advance <- function(prep, Y, t0, t1, h = 1e-4, rtol = 1e-6, atol = 1e-10) {
  u <- matrix(0, nrow(Y), ncol(Y))
  seg <- rosenbrock_segment(prep$f, Y, t0, t1, prep$V, u, h, rtol, atol,
                            uses_t = mrg_uses_time(prep$model))
  rownames(seg$Y) <- rownames(Y)
  seg
}
