# ============================================================================
# PBPK Simulator
#
# A whole-body, perfusion-limited PBPK model with mechanistic tissue
# partitioning (Rodgers & Rowland), hepatic and renal clearance, oral
# absorption with gut-wall extraction, and first-order pediatric scaling
# (allometric flows, GFR maturation, enzyme ontogeny).
#
# The structural model lives in models/pbpk.cpp (mrgsolve).
# ============================================================================

# Two engines, same answers. Run locally with mrgsolve installed, the app
# simulates through models/pbpk.cpp. In the browser (webR) nothing can be
# compiled, so it uses the exact matrix-exponential solution in
# R/pbpk_engine.R instead. tests/test_engine_vs_mrgsolve.R holds the two to
# within 1e-6 and the deploy fails if they drift apart.
library(shiny)
library(shinydashboard)
library(DT)
library(ggplot2)
library(scales)

# Locally this app runs from apps/small-molecule/ (shared code and models two
# levels up); in the browser build everything sits next to app.R.
source(if (file.exists(file.path("R", "theme.R"))) file.path("R", "theme.R") else file.path("..", "..", "shared", "theme.R"))
for (f in c("physiology", "partition", "drugs", "pbpk_engine", "simulate")) {
  source(file.path("R", paste0(f, ".R")))
}

ENGINE <- list(name = "Matrix exponential (base R)", simulate = pbpk_simulate)
REF <- file.path("..", "..", "reference", "pbpk_small_molecule.R")
if (file.exists(REF)) {
  source(REF)
  if (mrgsolve_ready()) {
    ok <- tryCatch({ pbpk_mrgsolve_model(file.path("..", "..", "models")); TRUE }, error = function(e) FALSE)
    if (ok) {
      ENGINE <- list(name = paste("mrgsolve", utils::packageVersion("mrgsolve")),
                     simulate = pbpk_simulate_mrgsolve)
    }
  }
}

plain_number <- function(x) format(x, scientific = FALSE, drop0trailing = TRUE, trim = TRUE)

muted <- function(...) tags$p(..., style = "color: #6b7078; font-size: 12px;")

TISSUE_LABELS <- c(lung = "Lung", brain = "Brain", heart = "Heart", kidney = "Kidney",
                   muscle = "Muscle", skin = "Skin", adipose = "Adipose", bone = "Bone",
                   spleen = "Spleen", gut = "Gut wall", liver = "Liver", rest = "Rest of body")

ROUTES <- c("Oral" = "oral", "IV bolus" = "iv_bolus", "IV infusion" = "iv_inf")

# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------

ui <- dashboardPage(
  skin = "blue",

  dashboardHeader(
    title = tags$div(
      tags$span("PBPK", style = "font-weight: bold; font-size: 24px;"),
      tags$span(" Simulator", style = "font-size: 18px; opacity: 0.9;")
    ),
    titleWidth = 350
  ),

  dashboardSidebar(
    width = 320,
    custom_css,
    sidebarMenu(
      menuItem("Plasma exposure", tabName = "dashboard", icon = icon("chart-line")),
      menuItem("Tissues", tabName = "tissues", icon = icon("lungs")),
      menuItem("Age scaling", tabName = "age", icon = icon("child")),
      menuItem("Drug & model setup", tabName = "setup", icon = icon("sliders")),
      menuItem("About", tabName = "about", icon = icon("info-circle"))
    ),

    tags$div(
      class = "sidebar-inputs",
      style = "padding: 18px; padding-top: 8px;",

      tags$h4("Compound", style = "font-size: 15px;"),
      selectInput("drug", NULL, choices = drug_choices, selected = "midazolam", width = "100%"),

      tags$h4("Subject", style = "font-size: 15px; margin-top: 10px;"),
      fluidRow(
        column(6, numericInput("age", "Age (years)", value = 30, min = 0, max = 90, step = 1)),
        column(6, numericInput("wt", "Weight (kg)", value = 70, min = 0.5, max = 200, step = 1))
      ),
      muted("Weight is refilled with the growth-chart median whenever the age changes; overwrite it if you like."),

      tags$h4("Dosing regimen", style = "font-size: 15px; margin-top: 10px;"),
      radioButtons("route", NULL, choices = ROUTES, selected = "oral", inline = TRUE),
      fluidRow(
        column(6, numericInput("dose", "Dose", value = 7.5, min = 0, step = 0.5)),
        column(6, selectInput("dose_basis", "Unit", choices = c("mg" = "flat", "mg/kg" = "mgkg"), selectize = FALSE))
      ),
      fluidRow(
        column(6, numericInput("interval", "Interval (h)", value = 24, min = 0.5, step = 1)),
        column(6, numericInput("n_doses", "Doses", value = 1, min = 1, max = 60, step = 1))
      ),
      conditionalPanel(
        "input.route == 'iv_inf'",
        numericInput("inf_dur", "Infusion duration (h)", value = 1, min = 0.05, step = 0.25, width = "100%")
      ),
      numericInput("duration", "Simulation duration (h)", value = 24, min = 1, step = 6, width = "100%")
    )
  ),

  dashboardBody(
    custom_css,
    tabItems(

      # --- Plasma exposure -------------------------------------------------
      tabItem(
        tabName = "dashboard",
        fluidRow(
          valueBoxOutput("cmax_box", width = 3),
          valueBoxOutput("auc_box", width = 3),
          valueBoxOutput("cl_box", width = 3),
          valueBoxOutput("thalf_box", width = 3)
        ),
        fluidRow(
          box(
            title = tags$div(
              tags$strong("Plasma concentration-time profile"),
              tags$span(textOutput("plot_subtitle", inline = TRUE),
                        style = "color: #6b7078; font-size: 14px; font-weight: normal;")
            ),
            status = "primary", solidHeader = TRUE, width = 12, collapsible = TRUE,
            checkboxInput("log_plasma", "Log concentration axis", value = FALSE),
            plotOutput("plasma_plot", height = "520px")
          )
        ),
        fluidRow(
          box(
            title = "Predicted disposition", status = "info", solidHeader = TRUE,
            width = 6, collapsible = TRUE,
            DT::dataTableOutput("disposition_table")
          ),
          box(
            title = "Exposure over the assessment window", status = "success",
            solidHeader = TRUE, width = 6, collapsible = TRUE,
            DT::dataTableOutput("exposure_table"),
            muted("The window runs from the last dose to one interval later, or to the end of the simulation for a single dose. Values are medians with the 5th-95th percentile across virtual subjects.")
          )
        )
      ),

      # --- Tissues -----------------------------------------------------------
      tabItem(
        tabName = "tissues",
        fluidRow(
          box(
            title = "Tissue concentrations, typical subject", status = "primary",
            solidHeader = TRUE, width = 12,
            checkboxInput("log_tissue", "Log concentration axis", value = TRUE),
            plotOutput("tissue_plot", height = "640px"),
            muted("Blue is the tissue; the grey dashed line is venous plasma for reference. Each panel has its own y-axis. Always computed with the matrix-exponential engine, which returns every compartment.")
          )
        ),
        fluidRow(
          box(
            title = "Tissue:plasma partition coefficients (Rodgers & Rowland)",
            status = "info", solidHeader = TRUE, width = 6,
            plotOutput("kp_plot", height = "380px")
          ),
          box(
            title = "Where the drug sits at steady state", status = "info",
            solidHeader = TRUE, width = 6,
            plotOutput("vd_plot", height = "380px"),
            muted("Each tissue's share of the volume of distribution, V_tissue x Kp. This is what drives the terminal half-life.")
          )
        )
      ),

      # --- Age scaling ---------------------------------------------------------
      tabItem(
        tabName = "age",
        fluidRow(
          box(
            title = "Clearance per kg across childhood", status = "primary",
            solidHeader = TRUE, width = 6,
            plotOutput("age_cl_plot", height = "400px"),
            muted("Typical subject at the growth-chart median weight for each age. Exact values from the rate matrix, not read off simulated curves.")
          ),
          box(
            title = "Profiles by age after the same mg/kg dose", status = "primary",
            solidHeader = TRUE, width = 6,
            plotOutput("age_profile_plot", height = "400px"),
            muted(textOutput("age_dose_note", inline = TRUE))
          )
        ),
        fluidRow(
          box(
            title = "Exposure-matched pediatric doses", status = "success",
            solidHeader = TRUE, width = 12,
            DT::dataTableOutput("age_table"),
            muted("Kinetics are linear, so the dose giving the adult AUC is the adult mg/kg dose times the ratio of AUC per mg/kg. The ontogeny curves are illustrative (see Drug & model setup); treat these as a demonstration of the method, not as dosing advice.")
          )
        )
      ),

      # --- Setup ---------------------------------------------------------------
      tabItem(
        tabName = "setup",
        fluidRow(
          box(
            title = "Physicochemistry and binding", status = "primary",
            solidHeader = TRUE, width = 4,
            selectInput("type", "Ionisation class",
                        choices = c("Neutral" = "neutral", "Monoprotic acid" = "acid",
                                    "Monoprotic base" = "base")),
            fluidRow(
              column(6, numericInput("logP", "log P", value = 3.13, step = 0.1)),
              column(6, numericInput("pKa", "pKa", value = 6, step = 0.1))
            ),
            fluidRow(
              column(6, numericInput("fu", "fu, plasma", value = 0.032, min = 0.0001, max = 1, step = 0.01)),
              column(6, numericInput("BP", "Blood:plasma", value = 0.664, min = 0.3, max = 5, step = 0.05))
            ),
            muted("Bases with pKa >= 7 bind acidic phospholipids, with the affinity back-calculated from B:P. Everything else binds albumin (acids, weak bases) or lipoproteins (neutrals).")
          ),
          box(
            title = "Clearance and absorption", status = "primary",
            solidHeader = TRUE, width = 4,
            fluidRow(
              column(6, numericInput("clint", "Hepatic CLint,u (L/h)", value = 1400, min = 0, step = 10)),
              column(6, numericInput("renal_factor", "Renal factor", value = 0, min = 0, max = 5, step = 0.05))
            ),
            selectInput("pathway", "Enzyme ontogeny", choices = names(ONTOGENY), selected = "CYP3A4"),
            fluidRow(
              column(4, numericInput("ka", "ka (1/h)", value = 3, min = 0.01, step = 0.1)),
              column(4, numericInput("fa", "Fa", value = 0.88, min = 0, max = 1, step = 0.05)),
              column(4, numericInput("fg", "Fg", value = 0.59, min = 0, max = 1, step = 0.05))
            ),
            muted("CLint is for a 70 kg adult and scales with liver mass and enzyme maturation. Renal clearance is fu x GFR x the renal factor (1 = filtration only; above 1 net secretion, below 1 net reabsorption).")
          ),
          box(
            title = "Variability and simulation", status = "info",
            solidHeader = TRUE, width = 4,
            fluidRow(
              column(4, numericInput("cv_clint", "CV CLint (%)", value = 35, min = 0, max = 150, step = 5)),
              column(4, numericInput("cv_co", "CV card. output (%)", value = 10, min = 0, max = 100, step = 5)),
              column(4, numericInput("cv_gfr", "CV GFR (%)", value = 15, min = 0, max = 100, step = 5))
            ),
            fluidRow(
              column(4, numericInput("n_subjects", "Subjects", value = 100, min = 1, max = 2000, step = 50)),
              column(4, numericInput("seed", "Seed", value = 123, min = 1, step = 1)),
              column(4, numericInput("delta", "Output step (h)", value = 0.1, min = 0.01, max = 2, step = 0.05))
            ),
            muted(textOutput("drug_note", inline = TRUE))
          )
        ),
        fluidRow(
          box(
            title = "Partition coefficients and tissue volumes, current subject",
            status = "info", solidHeader = TRUE, width = 6,
            DT::dataTableOutput("kp_table")
          ),
          box(
            title = "Enzyme and renal maturation", status = "info",
            solidHeader = TRUE, width = 6,
            plotOutput("ontogeny_plot", height = "360px"),
            muted("Fraction of adult activity per gram of liver (enzymes) and of adult GFR per 70 kg (Rhodin et al. 2009). Enzyme curves are illustrative Hill functions.")
          )
        )
      ),

      # --- About ---------------------------------------------------------------
      tabItem(
        tabName = "about",
        box(
          title = "About this application", status = "primary",
          solidHeader = TRUE, width = 12,
          tags$div(
            style = "padding: 20px;",
            tags$h3("PBPK Simulator"),
            tags$p("A whole-body physiologically based pharmacokinetic model you can
                    drive from a browser. Pick a compound or enter your own
                    physicochemical and clearance data, set the subject and the
                    regimen, and it predicts plasma and tissue concentrations,
                    exposure and how clearance changes across childhood."),
            tags$h4("Model structure"),
            tags$ul(
              tags$li("Twelve perfusion-limited tissues plus arterial and venous blood; spleen and gut drain into the liver through the portal vein"),
              tags$li("Tissue:plasma partition coefficients from the Rodgers & Rowland mechanistic equations"),
              tags$li("Hepatic metabolism of unbound drug in the liver at the intrinsic clearance; renal filtration of unbound drug"),
              tags$li("First-order oral absorption with fraction absorbed and gut-wall availability"),
              tags$li("Allometric cardiac output, GFR maturation (Rhodin 2009) and enzyme ontogeny for children"),
              tags$li("Log-normal variability on intrinsic clearance, cardiac output and GFR")
            ),
            tags$h4("Engine"),
            tags$p("This session is using: ", tags$strong(ENGINE$name), ". ",
                   "The model is defined in models/pbpk.cpp for mrgsolve. Because
                    it is linear, the browser build solves it exactly with a
                    matrix exponential instead; the test suite checks the two
                    against each other before every deployment."),
            tags$h4("Units"),
            tags$p("Amount mg, volume L, flow and clearance L/h, time h, concentration mg/L."),
            tags$hr(),
            tags$p(
              tags$strong("For research and teaching only. "),
              "The compound files are illustrations assembled from published
               inputs and calibrated clearances, and the pediatric ontogeny
               functions are simplified. Nothing here is validated for clinical
               use, and it must not be used to guide the treatment of a patient.",
              style = "color: #d1453b;"
            )
          )
        )
      )
    )
  )
)

# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------

server <- function(input, output, session) {

  # --- Presets --------------------------------------------------------------

  observeEvent(input$drug, {
    d <- DRUGS[[input$drug]]
    updateSelectInput(session, "type", selected = d$type)
    updateNumericInput(session, "logP", value = d$logP)
    updateNumericInput(session, "pKa", value = if (is.na(d$pKa)) 7 else d$pKa)
    updateNumericInput(session, "fu", value = d$fu)
    updateNumericInput(session, "BP", value = d$BP)
    updateNumericInput(session, "clint", value = d$clint)
    updateNumericInput(session, "renal_factor", value = d$renal_factor)
    updateSelectInput(session, "pathway", selected = d$pathway)
    updateNumericInput(session, "ka", value = d$ka)
    updateNumericInput(session, "fa", value = d$fa)
    updateNumericInput(session, "fg", value = d$fg)
    updateRadioButtons(session, "route", selected = d$route)
    updateNumericInput(session, "dose", value = d$dose)
    updateSelectInput(session, "dose_basis", selected = d$dose_basis)
    updateNumericInput(session, "interval", value = d$interval)
    updateNumericInput(session, "n_doses", value = d$n_doses)
    updateNumericInput(session, "duration", value = d$duration)
  })

  observeEvent(input$age, {
    shiny::req(is.finite(input$age))
    updateNumericInput(session, "wt", value = round(typical_weight(input$age), 1))
  }, ignoreInit = TRUE)

  output$drug_note <- renderText(DRUGS[[input$drug]]$note)

  # --- Inputs -> model objects ---------------------------------------------

  drug <- reactive({
    shiny::req(input$logP, input$fu, input$BP, input$clint, input$ka)
    validate(
      need(input$fu > 0 && input$fu <= 1, "fu must be between 0 and 1."),
      need(input$BP > 0, "Blood:plasma ratio must be positive."),
      need(input$clint >= 0, "CLint cannot be negative."),
      need(input$type == "neutral" || is.finite(input$pKa), "Enter a pKa for an ionisable compound.")
    )
    list(type = input$type, logP = input$logP,
         pKa = if (input$type == "neutral") NA else input$pKa,
         fu = input$fu, BP = input$BP, clint = input$clint,
         pathway = input$pathway, renal_factor = input$renal_factor,
         ka = input$ka, fa = input$fa, fg = input$fg)
  }) |> debounce(300)

  subject <- reactive({
    shiny::req(input$age, input$wt)
    validate(need(input$wt > 0, "Weight must be positive."),
             need(input$age >= 0, "Age cannot be negative."))
    list(age = input$age, wt = input$wt)
  }) |> debounce(300)

  regimen <- reactive({
    shiny::req(input$dose, input$interval, input$n_doses, input$duration)
    validate(
      need(input$dose > 0, "Dose must be positive."),
      need(input$interval > 0, "Dosing interval must be positive."),
      need(input$duration > 0, "Simulation duration must be positive."),
      need(input$route != "iv_inf" || (isTRUE(input$inf_dur > 0) && input$inf_dur <= input$interval),
           "Infusion duration must be positive and no longer than the dosing interval.")
    )
    amt <- if (input$dose_basis == "mgkg") input$dose * subject()$wt else input$dose
    list(route = input$route, amt = amt, interval = input$interval,
         n_doses = as.integer(input$n_doses),
         inf_dur = if (input$route == "iv_inf") input$inf_dur else 0)
  }) |> debounce(300)

  times <- reactive({
    shiny::req(input$delta, input$duration)
    validate(need(input$delta > 0, "Output step must be positive."))
    seq(0, input$duration, by = min(input$delta, input$duration / 20))
  })

  typical <- reactive({
    s <- subject()
    pbpk_parameters(drug(), s$wt, s$age)
  })

  metrics <- reactive(pbpk_metrics(typical()))

  population <- reactive({
    s <- subject()
    n <- input$n_subjects
    shiny::req(n)
    validate(need(n >= 1 && n <= 2000, "Use between 1 and 2000 virtual subjects."))
    make_population(drug(), s$wt, s$age, n, input$cv_clint, input$cv_co,
                    input$cv_gfr, input$seed)
  })

  sim <- reactive({
    reg <- regimen()
    tt <- times()
    res <- withProgress(message = "Simulating", value = 0.5,
                        ENGINE$simulate(population(), reg, tt))
    ex <- exposure_metrics(tt, res$plasma, reg)
    list(summary = summarise_profiles(tt, res$plasma), exposure = ex, reg = reg)
  })

  # --- Value boxes ------------------------------------------------------------

  stat_box <- function(value, unit, label, sub, icon_name, color = "blue") {
    valueBox(
      value = tags$div(
        tags$span(value, style = "font-size: 34px; font-weight: bold;"),
        tags$span(unit, style = "font-size: 18px;")
      ),
      subtitle = tags$div(tags$strong(label), tags$br(),
                          tags$span(sub, style = "font-size: 11px;")),
      icon = icon(icon_name), color = color, width = NULL
    )
  }

  fmt <- function(x) formatC(x, digits = 3, format = "fg", flag = "#")

  output$cmax_box <- renderValueBox({
    e <- sim()$exposure$per_subject
    stat_box(fmt(stats::median(e$cmax)), " mg/L", "Cmax", "Median, assessment window", "arrow-up")
  })

  output$auc_box <- renderValueBox({
    e <- sim()$exposure$per_subject
    stat_box(fmt(stats::median(e$auc)), " mg·h/L", "AUC", "Median, assessment window", "chart-area")
  })

  output$cl_box <- renderValueBox({
    m <- metrics()
    stat_box(fmt(m$CL), " L/h", "Plasma clearance",
             sprintf("Typical subject, %.3g L/h/kg", m$CL / subject()$wt), "tint")
  })

  output$thalf_box <- renderValueBox({
    m <- metrics()
    stat_box(fmt(m$t_half), " h", "Terminal half-life",
             sprintf("Vss %.3g L/kg", m$Vss / subject()$wt), "hourglass-half")
  })

  output$plot_subtitle <- renderText({
    sprintf(" (n = %s virtual subjects, %s)", format(input$n_subjects, big.mark = ","), ENGINE$name)
  })

  # --- Plasma plot -------------------------------------------------------------

  output$plasma_plot <- renderPlot({
    s <- sim()
    d <- s$summary
    win <- s$exposure$window
    log_y <- isTRUE(input$log_plasma)
    if (log_y) {
      # Pre-dose zeros have no place on a log axis; drop them rather than
      # clamping them to a floor, which would draw a spike at t = 0.
      floor_c <- max(d$med, na.rm = TRUE) * 1e-6
      d[, -1] <- lapply(d[, -1], function(v) ifelse(v > floor_c, v, NA))
    }
    # For a single dose the window is the whole simulation, so there's
    # nothing to mark.
    shade <- if (s$reg$n_doses > 1) {
      annotate("rect", xmin = win[1], xmax = win[2], ymin = -Inf, ymax = Inf,
               fill = PAL$sunken, colour = NA)
    }
    p <- ggplot(d, aes(x = time)) +
      shade +
      geom_ribbon(aes(ymin = q05, ymax = q95), fill = PAL$blue, alpha = 0.14, na.rm = TRUE) +
      geom_ribbon(aes(ymin = q25, ymax = q75), fill = PAL$blue, alpha = 0.28, na.rm = TRUE) +
      geom_line(aes(y = med), colour = PAL$blue_ink, linewidth = 1.1, na.rm = TRUE) +
      labs(x = "Time (h)", y = "Plasma concentration (mg/L)",
           title = sprintf("%s, %s", DRUGS[[input$drug]]$label, names(ROUTES)[ROUTES == input$route]),
           subtitle = paste0("Median with 50% and 90% prediction intervals.",
                             if (s$reg$n_doses > 1) " The shaded band is the assessment window." else "")) +
      scale_x_continuous(breaks = scales::pretty_breaks(10), expand = expansion(c(0, 0.02))) +
      theme_sim()
    if (log_y) p + scale_y_log10(labels = plain_number) else p + scale_y_continuous(limits = c(0, NA), expand = expansion(c(0, 0.05)))
  })

  # --- Tables ------------------------------------------------------------------

  small_table <- function(df) {
    DT::datatable(df, options = list(dom = "t", ordering = FALSE, pageLength = 50),
                  rownames = FALSE) |>
      DT::formatStyle(1, fontWeight = "bold", color = "#16181d")
  }

  output$disposition_table <- DT::renderDataTable({
    m <- metrics(); p <- typical(); s <- subject()
    small_table(data.frame(
      Parameter = c("Plasma clearance", "Hepatic extraction ratio", "Oral bioavailability (Fa x Fg x Fh)",
                    "Volume at steady state", "Terminal half-life", "Hepatic CLint,u (this subject)",
                    "Enzyme activity vs adult", "Renal clearance", "GFR", "Cardiac output"),
      Value = c(sprintf("%.3g L/h (%.3g L/h/kg)", m$CL, m$CL / s$wt),
                sprintf("%.2f", m$E_H),
                sprintf("%.2f", m$F_oral),
                sprintf("%.3g L (%.3g L/kg)", m$Vss, m$Vss / s$wt),
                sprintf("%.3g h", m$t_half),
                sprintf("%.4g L/h", m$CLint),
                sprintf("%.0f%%", 100 * enzyme_maturation(s$age, input$pathway)),
                sprintf("%.3g L/h", m$CLr),
                sprintf("%.3g L/h (%.0f%% of adult per 70 kg)", m$GFR, 100 * gfr_maturation(s$age)),
                sprintf("%.3g L/h", p$CO))
    ))
  })

  output$exposure_table <- DT::renderDataTable({
    e <- sim()$exposure$per_subject
    row <- function(x, unit) sprintf("%.3g (%.3g - %.3g) %s", stats::median(x),
                                     stats::quantile(x, 0.05), stats::quantile(x, 0.95), unit)
    small_table(data.frame(
      Metric = c("Cmax", "Time of Cmax after last dose", "Concentration at end of window", "AUC over window"),
      Value = c(row(e$cmax, "mg/L"), row(e$tmax, "h"), row(e$cmin, "mg/L"), row(e$auc, "mg·h/L"))
    ))
  })

  # --- Tissues -----------------------------------------------------------------

  tissue_sim <- reactive({
    tt <- times()
    pbpk_simulate(list(typical()), regimen(), tt, keep_states = TRUE)
  })

  output$tissue_plot <- renderPlot({
    r <- tissue_sim()
    st <- r$states[, , 1]
    long <- do.call(rbind, lapply(names(TISSUE_LABELS), function(t) {
      data.frame(time = r$time, tissue = TISSUE_LABELS[[t]], conc = st[, toupper(t)])
    }))
    long$tissue <- factor(long$tissue, levels = TISSUE_LABELS)
    plasma <- data.frame(time = r$time, conc = r$plasma[, 1])
    log_y <- isTRUE(input$log_tissue)
    if (log_y) {
      floor_c <- max(long$conc) * 1e-6
      long$conc <- ifelse(long$conc > floor_c, long$conc, NA)
      plasma$conc <- ifelse(plasma$conc > floor_c, plasma$conc, NA)
    }
    p <- ggplot(long, aes(time, conc)) +
      geom_line(data = plasma, colour = PAL$ink_3, linetype = "22", linewidth = 0.6, na.rm = TRUE) +
      geom_line(colour = PAL$blue, linewidth = 0.9, na.rm = TRUE) +
      facet_wrap(~tissue, ncol = 4, scales = "free_y") +
      labs(x = "Time (h)", y = "Concentration (mg/L)",
           title = "Tissue concentrations",
           subtitle = sprintf("Typical %s-year-old, %s kg. Dashed: venous plasma.", format(subject()$age), format(subject()$wt))) +
      theme_sim(12) +
      theme(strip.text = element_text(colour = PAL$ink, face = "bold", hjust = 0),
            panel.spacing = unit(14, "pt"))
    if (log_y) p + scale_y_log10(labels = plain_number) else p
  })

  output$kp_plot <- renderPlot({
    p <- typical()
    d <- data.frame(tissue = TISSUE_LABELS[names(p$Kp)], kp = unname(p$Kp))
    d$tissue <- stats::reorder(d$tissue, d$kp)
    ggplot(d, aes(kp, tissue)) +
      geom_col(fill = PAL$blue, width = 0.7) +
      geom_text(aes(label = sprintf("%.2f", kp)), hjust = -0.15, colour = PAL$ink_2, size = 3.6) +
      geom_vline(xintercept = 1, colour = PAL$ink_3, linetype = "22") +
      scale_x_continuous(expand = expansion(c(0, 0.15))) +
      labs(x = "Kp (tissue:plasma)", y = NULL, title = NULL) +
      theme_sim(12) + theme(panel.grid.major.y = element_blank())
  })

  output$vd_plot <- renderPlot({
    p <- typical()
    vt <- p$V[names(p$Kp)] * p$Kp
    blood <- (p$V[["art"]] + p$V[["ven"]]) * p$BP
    d <- data.frame(tissue = c(TISSUE_LABELS[names(p$Kp)], "Blood"), v = c(unname(vt), blood))
    d$share <- d$v / sum(d$v)
    d$tissue <- stats::reorder(d$tissue, d$share)
    ggplot(d, aes(share, tissue)) +
      geom_col(fill = PAL$blue, width = 0.7) +
      geom_text(aes(label = sprintf("%.0f%%", 100 * share)), hjust = -0.15, colour = PAL$ink_2, size = 3.6) +
      scale_x_continuous(labels = scales::percent, expand = expansion(c(0, 0.15))) +
      labs(x = "Share of the volume of distribution", y = NULL) +
      theme_sim(12) + theme(panel.grid.major.y = element_blank())
  })

  # --- Age scaling --------------------------------------------------------------

  ref_dose_mgkg <- reactive({
    if (input$dose_basis == "mgkg") input$dose else input$dose / REF_WT
  })

  output$age_cl_plot <- renderPlot({
    d <- drug()
    ages <- c(seq(1 / 365, 1, length.out = 40), seq(1.1, 18, length.out = 40))
    tab <- age_scaling_table(d, ages, 1, oral = FALSE)
    key <- age_scaling_table(d, AGE_LABELS[AGE_LABELS < 30], 1, oral = FALSE)$table
    ggplot(tab$table, aes(age, CL_kg)) +
      geom_hline(yintercept = tab$adult$CL_kg, colour = PAL$ink_3, linetype = "22") +
      annotate("text", x = 1 / 365, y = tab$adult$CL_kg, label = "Adult (70 kg)", vjust = -0.6, hjust = 0,
               colour = PAL$ink_3, size = 3.6) +
      geom_line(colour = PAL$blue, linewidth = 1.1) +
      geom_point(data = key, colour = PAL$blue_ink, size = 2.6) +
      scale_x_log10(breaks = c(1 / 52, 1 / 12, 0.25, 1, 2, 6, 12, 18),
                    labels = c("1 wk", "1 mo", "3 mo", "1 y", "2 y", "6 y", "12 y", "18 y")) +
      scale_y_continuous(limits = c(0, NA), expand = expansion(c(0, 0.08))) +
      labs(x = "Postnatal age (log scale)", y = "Plasma clearance (L/h/kg)",
           title = NULL) +
      theme_sim(12)
  })

  age_profiles <- reactive({
    d <- drug()
    ages <- AGE_LABELS[c("1 week", "3 months", "1 year", "6 years", "12 years", "Adult")]
    tt <- times()
    reg <- regimen()
    do.call(rbind, lapply(seq_along(ages), function(i) {
      a <- ages[[i]]
      w <- if (names(ages)[i] == "Adult") REF_WT else typical_weight(a)
      reg$amt <- ref_dose_mgkg() * w
      r <- pbpk_simulate(list(pbpk_parameters(d, w, a)), reg, tt)
      data.frame(time = tt, conc = r$plasma[, 1], age = names(ages)[i])
    }))
  })

  output$age_dose_note <- renderText({
    sprintf("%.3g mg/kg for every age (the sidebar dose%s), same route and schedule.",
            ref_dose_mgkg(), if (input$dose_basis == "flat") " divided by 70 kg" else "")
  })

  output$age_profile_plot <- renderPlot({
    d <- age_profiles()
    lv <- unique(d$age)
    d$age <- factor(d$age, levels = lv)
    # Log axis: the ages differ mostly in how fast the drug is eliminated,
    # which is a difference in slope.
    d$conc <- ifelse(d$conc > max(d$conc) * 1e-6, d$conc, NA)
    ggplot(d, aes(time, conc, colour = age)) +
      geom_line(linewidth = 1, na.rm = TRUE) +
      scale_colour_manual(values = stats::setNames(SERIES[seq_along(lv)], lv)) +
      scale_y_log10(labels = plain_number) +
      labs(x = "Time (h)", y = "Plasma concentration (mg/L, log scale)") +
      theme_sim(12)
  })

  output$age_table <- DT::renderDataTable({
    oral <- input$route == "oral"
    res <- age_scaling_table(drug(), AGE_LABELS, ref_dose_mgkg(), oral)
    t <- res$table
    df <- data.frame(
      Age = names(AGE_LABELS),
      `Weight (kg)` = sprintf("%.1f", t$wt),
      `Enzyme vs adult` = sprintf("%.0f%%", 100 * t$enzyme),
      `GFR maturation` = sprintf("%.0f%%", 100 * t$gfr_mat),
      `CL (L/h/kg)` = sprintf("%.3g", t$CL_kg),
      `Vss (L/kg)` = sprintf("%.3g", t$Vss_kg),
      `t½ (h)` = sprintf("%.3g", t$t_half),
      `AUC at this dose (mg·h/L)` = sprintf("%.3g", t$auc),
      `Dose for adult AUC (mg/kg)` = sprintf("%.3g", t$dose_match),
      check.names = FALSE
    )
    DT::datatable(df, options = list(dom = "t", ordering = FALSE, pageLength = 20), rownames = FALSE) |>
      DT::formatStyle(1, fontWeight = "bold", color = "#16181d")
  })

  # --- Setup outputs -------------------------------------------------------------

  output$kp_table <- DT::renderDataTable({
    p <- typical()
    tis <- names(p$Kp)
    small_table(data.frame(
      Tissue = TISSUE_LABELS[tis],
      `Volume (L)` = sprintf("%.3g", p$V[tis]),
      `Blood flow (L/h)` = sprintf("%.3g", c(p$Q[setdiff(tis, "lung")], lung = p$CO)[tis]),
      Kp = sprintf("%.3g", p$Kp),
      check.names = FALSE
    ))
  })

  output$ontogeny_plot <- renderPlot({
    ages <- 10^seq(log10(1 / 365), log10(18), length.out = 120)
    d <- rbind(
      do.call(rbind, lapply(setdiff(names(ONTOGENY), "Mature (no ontogeny)"), function(k) {
        data.frame(age = ages, f = enzyme_maturation(ages, k), what = k)
      })),
      data.frame(age = ages, f = gfr_maturation(ages), what = "GFR")
    )
    lv <- c("CYP3A4", "CYP1A2", "CYP2D6", "Generic hepatic", "GFR")
    d$what <- factor(d$what, levels = lv)
    ggplot(d, aes(age, f, colour = what)) +
      geom_line(aes(linewidth = what == input$pathway)) +
      scale_linewidth_manual(values = c(`TRUE` = 1.6, `FALSE` = 0.8), guide = "none") +
      scale_colour_manual(values = stats::setNames(SERIES[seq_along(lv)], lv)) +
      scale_x_log10(breaks = c(1 / 52, 1 / 12, 0.25, 1, 2, 6, 18),
                    labels = c("1 wk", "1 mo", "3 mo", "1 y", "2 y", "6 y", "18 y")) +
      scale_y_continuous(labels = scales::percent, limits = c(0, 1.05)) +
      labs(x = "Postnatal age (log scale)", y = "Fraction of adult") +
      theme_sim(12)
  })
}

shinyApp(ui, server)
