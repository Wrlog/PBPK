# ============================================================================
# Appearance: one palette, shared by the page chrome and the plots.
#
# The colours are a colour-vision-safe set. Blue carries the simulated
# concentration, amber and red the two targets, and grey everything that is
# context rather than data. Nothing is a gradient and nothing animates: the
# panel should read like a result, not a product page.
# ============================================================================

PAL <- list(
  ink      = "#16181d",
  ink_2    = "#454951",
  ink_3    = "#6b7078",
  rule     = "#e4e6ea",
  surface  = "#ffffff",
  sunken   = "#f7f8fa",
  blue     = "#2a78d6",
  blue_ink = "#1d5fae",
  blue_bg  = "#edf3fc",
  green    = "#1baf7a",
  green_ink= "#0d7a54",
  amber    = "#eda100",
  amber_ink= "#8a6200",
  red      = "#d1453b",
  violet   = "#4a3aa7"
)

app_css <- sprintf("
  :root {
    --ink: %s; --ink-2: %s; --ink-3: %s; --rule: %s;
    --surface: %s; --sunken: %s;
    --blue: %s; --blue-ink: %s; --blue-bg: %s;
  }

  body, .content-wrapper, .right-side {
    background: var(--sunken);
    font-family: ui-sans-serif, -apple-system, 'Segoe UI', Roboto, sans-serif;
    color: var(--ink-2);
    font-size: 15px;
  }

  /* --- header ----------------------------------------------------------- */

  .skin-blue .main-header .navbar,
  .skin-blue .main-header .logo,
  .skin-blue .main-header .logo:hover {
    background: var(--surface);
    border-bottom: 1px solid var(--rule);
    box-shadow: none;
  }
  .skin-blue .main-header .logo {
    color: var(--ink);
    font-weight: 650;
    font-size: 17px;
    letter-spacing: -0.01em;
    border-right: 1px solid var(--rule);
  }
  .skin-blue .main-header .navbar .sidebar-toggle { color: var(--ink-3); }
  .skin-blue .main-header .navbar .sidebar-toggle:hover {
    background: var(--sunken); color: var(--ink);
  }

  /* --- sidebar: light, not a dark slab ---------------------------------- */

  .skin-blue .main-sidebar, .skin-blue .wrapper {
    background: var(--surface);
  }
  .skin-blue .main-sidebar { border-right: 1px solid var(--rule); }
  .sidebar { padding-bottom: 40px; }
  .sidebar .shiny-input-container { margin-bottom: 16px; color: var(--ink-2); }
  .sidebar .shiny-input-container > label,
  .sidebar .control-label {
    color: var(--ink-2) !important;
    font-weight: 500;
    font-size: 13px;
    margin-bottom: 5px;
  }
  .sidebar h3, .sidebar h4 {
    color: var(--ink-3);
    font-size: 11px;
    font-weight: 700;
    text-transform: uppercase;
    letter-spacing: 0.09em;
    margin: 22px 0 12px;
    padding-bottom: 8px;
    border-bottom: 1px solid var(--rule);
  }
  .sidebar input[type=number], .sidebar input[type=text], .sidebar select {
    background: var(--surface);
    color: var(--ink);
    border: 1px solid var(--rule);
    border-radius: 7px;
    padding: 6px 9px;
    width: 100%%;
    box-shadow: none;
  }
  .sidebar input[type=number]:focus, .sidebar select:focus {
    border-color: var(--blue);
    outline: 2px solid rgba(42,120,214,0.18);
    outline-offset: 0;
  }
  .sidebar .radio label, .sidebar .checkbox label,
  .sidebar .shiny-options-group label {
    color: var(--ink-2) !important;
    font-weight: 400;
    font-size: 13.5px;
  }

  /* --- sliders ---------------------------------------------------------- */

  .irs--shiny .irs-bar { background: var(--blue); border: none; }
  .irs--shiny .irs-line { background: var(--rule); border: none; }
  .irs--shiny .irs-handle {
    border: 2px solid var(--blue);
    background: var(--surface);
    box-shadow: none;
    width: 18px; height: 18px; top: 22px;
  }
  .irs--shiny .irs-single, .irs--shiny .irs-from, .irs--shiny .irs-to {
    background: var(--blue-ink); border-radius: 5px; font-size: 11px;
  }
  .irs--shiny .irs-min, .irs--shiny .irs-max {
    background: transparent; color: var(--ink-3); font-size: 11px;
  }

  /* --- boxes ------------------------------------------------------------ */

  .box {
    border: 1px solid var(--rule);
    border-top: 1px solid var(--rule);
    border-radius: 10px;
    box-shadow: none;
    background: var(--surface);
  }
  .box-header { border-bottom: 1px solid var(--rule); padding: 14px 18px; }
  .box-header .box-title {
    font-size: 14px; font-weight: 650; color: var(--ink);
    letter-spacing: -0.005em;
  }
  .box-body { padding: 18px; }
  .box.box-solid > .box-header, .box.box-primary > .box-header {
    background: var(--surface); color: var(--ink);
  }
  .box.box-solid, .box.box-primary { border-top-color: var(--rule); }

  /* --- value boxes: flat, quiet, legible -------------------------------- */

  .small-box {
    border-radius: 10px;
    border: 1px solid var(--rule);
    box-shadow: none;
    background: var(--surface) !important;
    color: var(--ink) !important;
  }
  .small-box .icon { display: none; }
  .small-box p, .small-box h3 { color: var(--ink) !important; }
  .small-box > .inner { padding: 16px 18px; }
  .small-box.bg-green { border-top: 3px solid %s; }
  .small-box.bg-yellow { border-top: 3px solid %s; }
  .small-box.bg-red { border-top: 3px solid %s; }
  .small-box.bg-blue, .small-box.bg-aqua { border-top: 3px solid %s; }

  /* --- tables ----------------------------------------------------------- */

  table.dataTable thead th {
    border-bottom: 1px solid var(--rule) !important;
    color: var(--ink-3);
    font-size: 11.5px;
    text-transform: uppercase;
    letter-spacing: 0.04em;
    font-weight: 600;
  }
  table.dataTable tbody td {
    border-top: 1px solid var(--rule);
    color: var(--ink-2);
    font-variant-numeric: tabular-nums;
    font-size: 13.5px;
  }
  .dataTables_wrapper .dataTables_info,
  .dataTables_wrapper .dataTables_paginate { color: var(--ink-3) !important; font-size: 12.5px; }

  /* --- misc ------------------------------------------------------------- */

  .shiny-notification { border-radius: 8px; border: 1px solid var(--rule); }
  .help-block, .shiny-output-error-validation {
    color: var(--ink-3); font-size: 13px;
  }
  hr { border-top: 1px solid var(--rule); }
",
  PAL$ink, PAL$ink_2, PAL$ink_3, PAL$rule,
  PAL$surface, PAL$sunken,
  PAL$blue, PAL$blue_ink, PAL$blue_bg,
  PAL$green, PAL$amber, PAL$red, PAL$blue
)

custom_css <- tags$head(
  tags$meta(name = "viewport", content = "width=device-width, initial-scale=1"),
  tags$style(HTML(app_css))
)

#' Shared ggplot theme
#'
#' Recessive grid, no panel border, labels in ink rather than in the series
#' colour, so the only saturated thing on the panel is the data.
theme_sim <- function(base_size = 13) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = base_size + 3, face = "bold",
                                         colour = PAL$ink,
                                         margin = ggplot2::margin(b = 4)),
      plot.subtitle = ggplot2::element_text(size = base_size - 1,
                                            colour = PAL$ink_3,
                                            margin = ggplot2::margin(b = 14)),
      axis.title = ggplot2::element_text(size = base_size - 1, colour = PAL$ink_3),
      axis.text = ggplot2::element_text(size = base_size - 2, colour = PAL$ink_3),
      panel.grid.major = ggplot2::element_line(colour = PAL$rule, linewidth = 0.4),
      panel.grid.minor = ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(fill = PAL$surface, colour = NA),
      panel.background = ggplot2::element_rect(fill = PAL$surface, colour = NA),
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      legend.position = "top",
      legend.justification = "left",
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = base_size - 2, colour = PAL$ink_2),
      plot.margin = ggplot2::margin(14, 16, 10, 10)
    )
}

# Categorical order for plots with several series (ages, dose levels).
# Assigned in this order and never cycled: a series keeps its colour when
# others are added or removed. Validated for colour-vision deficiency on
# adjacent pairs.
SERIES <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100",
            "#e87ba4", "#008300", "#4a3aa7", "#e34948")
