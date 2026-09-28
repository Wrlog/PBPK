# ============================================================================
# Pieces shared by every app in this repository: finding files whether the
# app runs from the repository or from its browser build, picking the
# simulation engine, and a few UI and plotting helpers.
# ============================================================================

# Locally an app runs from apps/<name>/ with shared code in ../../shared and
# models in ../../models. In the browser build everything is copied next to
# app.R (R/ and models/). repo_file() finds a file in either layout.
repo_file <- function(...) {
  local <- file.path(...)
  if (file.exists(local)) return(local)
  up <- file.path("..", "..", ...)
  if (file.exists(up)) return(up)
  local
}

#' Load a model and choose how to solve it
#'
#' Run locally with mrgsolve installed, simulations go through mrgsolve (the
#' reference). In the browser the reference/ directory isn't shipped and the
#' base-R engine in ode_engine.R solves the same model file.
load_model <- function(cpp) {
  model <- mrg_read(repo_file("models", cpp))
  engine <- list(name = "Base-R Rosenbrock (browser)", solve = mrg_solve, mrgsolve = FALSE)
  ref <- repo_file("reference", "mrgsolve_engine.R")
  if (file.exists(ref)) {
    source(ref)
    if (mrgsolve_ready()) {
      ok <- tryCatch({ mrgsolve_model(model$path); TRUE }, error = function(e) FALSE)
      if (ok) {
        engine <- list(name = paste("mrgsolve", utils::packageVersion("mrgsolve")),
                       solve = mrg_solve_mrgsolve, mrgsolve = TRUE)
      }
    }
  }
  list(model = model, engine = engine)
}

plain_number <- function(x) format(x, scientific = FALSE, drop0trailing = TRUE, trim = TRUE)

muted <- function(...) tags$p(..., style = "color: #6b7078; font-size: 12px;")

stat_box <- function(value, unit, label, sub, icon_name, color = "blue") {
  valueBox(
    value = tags$div(
      tags$span(value, style = "font-size: 32px; font-weight: bold;"),
      tags$span(unit, style = "font-size: 17px;")
    ),
    subtitle = tags$div(tags$strong(label), tags$br(),
                        tags$span(sub, style = "font-size: 11px;")),
    icon = icon(icon_name), color = color, width = NULL
  )
}

#' Median and prediction-interval bands across columns (subjects)
summarise_bands <- function(time, x) {
  q <- apply(x, 1, stats::quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
             na.rm = TRUE, names = FALSE)
  if (is.null(dim(q))) q <- matrix(q, nrow = 5)
  data.frame(time = time, q05 = q[1, ], q25 = q[2, ], med = q[3, ], q75 = q[4, ], q95 = q[5, ])
}

band_plot <- function(d, x_lab, y_lab, subtitle = NULL) {
  ggplot(d, aes(x = time)) +
    geom_ribbon(aes(ymin = q05, ymax = q95), fill = PAL$blue, alpha = 0.14, na.rm = TRUE) +
    geom_ribbon(aes(ymin = q25, ymax = q75), fill = PAL$blue, alpha = 0.28, na.rm = TRUE) +
    geom_line(aes(y = med), colour = PAL$blue_ink, linewidth = 1.1, na.rm = TRUE) +
    labs(x = x_lab, y = y_lab, subtitle = subtitle) +
    theme_sim(12)
}

#' Lines for several named series, coloured in the fixed categorical order
series_plot <- function(d, x_lab, y_lab, levels = unique(d$series), subtitle = NULL) {
  d$series <- factor(d$series, levels = levels)
  ggplot(d, aes(time, value, colour = series)) +
    geom_line(linewidth = 1, na.rm = TRUE) +
    scale_colour_manual(values = stats::setNames(SERIES[seq_along(levels)], levels)) +
    labs(x = x_lab, y = y_lab, subtitle = subtitle) +
    theme_sim(12)
}

small_table <- function(df, page = 50) {
  DT::datatable(df, options = list(dom = "t", ordering = FALSE, pageLength = page),
                rownames = FALSE) |>
    DT::formatStyle(1, fontWeight = "bold", color = "#16181d")
}

#' The About tab: what the model is, where it comes from, how it runs
about_tab <- function(title, intro, structure, reference, engine_name, extra = NULL) {
  tabItem(
    tabName = "about",
    box(
      title = "About this model", status = "primary", solidHeader = TRUE, width = 12,
      tags$div(
        style = "padding: 20px;",
        tags$h3(title),
        tags$p(intro),
        tags$h4("Model structure"),
        tags$ul(lapply(structure, tags$li)),
        tags$h4("Source"),
        tags$p(reference),
        extra,
        tags$h4("Engine"),
        tags$p("This session is using: ", tags$strong(engine_name), ". The model is an
                mrgsolve model file. Run locally with mrgsolve installed, the app
                simulates through mrgsolve. In the browser, where C++ can't be
                compiled, the same file is translated to R and solved with a stiff
                Rosenbrock solver; the test suite checks the two against each other."),
        tags$hr(),
        tags$p(tags$strong("For research and teaching only. "),
               "Nothing here is validated for clinical use, and it must not be used to
                guide the treatment of a patient.", style = "color: #d1453b;")
      )
    )
  )
}
