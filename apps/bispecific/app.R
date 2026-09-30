# ============================================================================
# Bispecific antibodies: a generalized minimal PBPK model with two targets
# (Spinosa, Joslyn, Ramanujan, Gadkar & Hosseini, CPT:PSP 2026)
#
# The minimal PBPK body (blood, leaky and tight tissue, lymph, and a tumour)
# with a drug that has two different arms. Each arm binds its own target, and
# a bound drug can pick up the second target on another cell (a T-cell
# engager), on the same cell (a tumour-targeted bispecific, where avidity
# matters), or neutralise two soluble cytokines. The four case studies of the
# paper are the presets.
# ============================================================================

library(shiny)
library(shinydashboard)
library(DT)
library(ggplot2)
library(scales)

shared <- function(f) source(if (file.exists(file.path("R", f))) file.path("R", f) else file.path("..", "..", "shared", f))
for (f in c("theme.R", "ode_engine.R", "app_helpers.R")) shared(f)
source(file.path("R", "cases.R"))

M <- load_model("bispecific_mpbpk_spinosa2026.cpp")

case_choices <- stats::setNames(names(CASES), vapply(CASES, `[[`, "", "label"))

ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = tags$div(tags$span("PBPK", style = "font-weight: bold; font-size: 24px;"),
                                   tags$span(" Bispecific antibodies", style = "font-size: 18px;")),
                  titleWidth = 350),
  dashboardSidebar(
    width = 320, custom_css,
    sidebarMenu(
      menuItem("Drug and targets", tabName = "main", icon = icon("vial")),
      menuItem("Dose response", tabName = "dose", icon = icon("chart-line")),
      menuItem("Model setup", tabName = "setup", icon = icon("sliders")),
      menuItem("About", tabName = "about", icon = icon("info-circle"))
    ),
    tags$div(
      style = "padding: 18px; padding-top: 8px;",
      tags$h4("Case study"),
      selectInput("case", NULL, choices = case_choices, selected = "tce", width = "100%"),
      uiOutput("case_controls"),
      tags$h4("Regimen"),
      uiOutput("regimen_controls"),
      numericInput("days", "Simulate (days)", 28, min = 7, max = 365, step = 7, width = "100%")
    )
  ),
  dashboardBody(
    custom_css,
    tabItems(
      tabItem(
        tabName = "main",
        fluidRow(valueBoxOutput("vb1", width = 3), valueBoxOutput("vb2", width = 3),
                 valueBoxOutput("vb3", width = 3), valueBoxOutput("vb4", width = 3)),
        fluidRow(
          box(title = tags$div(tags$strong("Drug concentration"),
                               tags$span(textOutput("engine_note", inline = TRUE),
                                         style = "color: #6b7078; font-size: 13px; font-weight: normal;")),
              status = "primary", solidHeader = TRUE, width = 6, plotOutput("pk_plot", height = "380px"),
              uiOutput("pk_note")),
          box(title = "Target engagement", status = "primary", solidHeader = TRUE, width = 6,
              plotOutput("te_plot", height = "380px"), uiOutput("te_note"))
        ),
        fluidRow(box(title = "About this case study", status = "primary", solidHeader = TRUE, width = 12,
                     uiOutput("case_note")))
      ),
      tabItem(
        tabName = "dose",
        box(title = textOutput("dose_title", inline = TRUE), status = "primary", solidHeader = TRUE, width = 12,
            plotOutput("dose_plot", height = "440px"), uiOutput("dose_note"))
      ),
      tabItem(
        tabName = "setup",
        fluidRow(
          box(title = "Drug", status = "primary", solidHeader = TRUE, width = 6,
              fluidRow(column(6, numericInput("CLp", "Nonspecific clearance (mL/day/kg)", 4, min = 0, step = 0.1)),
                       column(6, numericInput("MWab", "Molecular weight (g/mol)", 150000, min = 1e4, step = 1e4))),
              fluidRow(column(6, numericInput("Kd_R1", "KD, arm 1 (nM)", 68, min = 1e-4, step = 1)),
                       column(6, numericInput("Kd_R2", "KD, arm 2 (nM)", 40, min = 1e-4, step = 1))),
              fluidRow(column(6, numericInput("kon", "kon, both arms (1/nM/day)", 3.6, min = 0, step = 0.1)),
                       column(6, numericInput("BW", "Body weight (kg)", 3.12, min = 0.02, step = 1))),
              muted("Changing the case study resets these to its published values. KD = koff / kon.")),
          box(title = "All parameters of this case", status = "primary", solidHeader = TRUE, width = 6,
              DTOutput("par_table"),
              muted("From the paper's Model parameters.xlsx (Data S1), with the app's choices applied."))
        )
      ),
      about_tab(
        "Generalized minimal PBPK model for bispecific antibodies",
        "A minimal PBPK model - blood, leaky and tight tissue, lymph, and optional efficacy and
         safety tissues - in which a drug with two different arms binds two targets in every
         space except lymph. The targets can be soluble or on cells; a drug bound by one arm can
         bind the other target on a different cell (trans, as a T-cell engager does), on the same
         cell (cis, with an avidity factor), or, for a monospecific antibody, bind the same target
         twice. The same model file runs all four case studies of the paper.",
        list("Transport: convection from blood into leaky and tight tissue (and the tumour) through vascular pores with reflection coefficients; return through lymph (Cao & Jusko 2014)",
             "Clearance: nonspecific, from blood only; drug-soluble target complexes clear like free drug and, optionally, move with it",
             "Targets: zero-order synthesis and first-order degradation in every space, set so the total target stays at its baseline; complexes internalise at the target's degradation rate",
             "Second binding: rate multiplied by an effective avidity chi / (Vchi x cells per mL); for a trimer, internalisation runs through both targets and the other target is recycled",
             "Built from the authors' SimBiology project by tools/build_bispecific.py: 93 states, 273 reactions, 12 rate rules"),
        tags$span("Spinosa P, Joslyn L, Ramanujan S, Gadkar K, Hosseini I. A Generalized Minimal PBPK-PD
                   Model of Bispecific Antibodies: Case Studies and Applications in Drug Development.
                   CPT Pharmacometrics Syst Pharmacol 2026;15:e70167. Minimal PBPK: Cao Y, Jusko WJ.
                   J Pharmacokinet Pharmacodyn 2014;41:375-387."),
        M$engine$name,
        extra = tagList(
          tags$h4("Further reading"),
          tags$ul(
            tags$li("Schropp J, Khot A, Shah DK, Koch G. Target-mediated drug disposition model for bispecific antibodies: properties, approximation, and optimal dosing strategy. CPT:PSP 2019;8:177-187."),
            tags$li("Susilo ME, Schaller S, Jimenez-Franco LD, et al. Whole-body physiologically based pharmacokinetic modeling framework for tissue target engagement of CD3 bispecific antibodies. Pharmaceutics 2025;17:500."),
            tags$li("Betts A, van der Graaf PH. Mechanistic quantitative pharmacology strategies for the early clinical development of bispecific antibodies in oncology. Clin Pharmacol Ther 2020;108:528-541."),
            tags$li("Gibbs JP, Yuraszeck T, Biesdorf C, Xu Y, Kasichayanula S. Informing development of bispecific antibodies using physiologically based pharmacokinetic-pharmacodynamic models. J Clin Pharmacol 2020;60:S132-S146.")
          )
        )
      )
    )
  )
)

server <- function(input, output, session) {

  cs <- reactive(CASES[[input$case]])

  output$case_controls <- renderUI({
    switch(input$case,
      affinity = radioButtons("cd3", "CD3 arm", c("Moderate affinity (40 nM)" = "40", "Low affinity (400 nM)" = "400",
                                                  "No CD3 binding" = "none"), selected = "40"),
      cis = tagList(
        radioButtons("format", "Molecule", c("Tumour-targeted bispecific" = "bispecific",
                                             "Bivalent antibody (10 nM, therapeutic receptor only)" = "bivalent")),
        selectInput("chi", "Effective avidity of the second arm", c(1, 10, 100, 1000, 1e4, 1e5, 1e6), selected = 1000),
        muted("Avidity multiplies the on-rate of the second binding, once the drug is held on a cell by its first arm.")
      ),
      NULL)
  })

  output$regimen_controls <- renderUI({
    if (input$case == "soluble") {
      tagList(
        radioButtons("route", NULL, c("Subcutaneous" = "SC", "Intravenous" = "IV"), inline = TRUE),
        fluidRow(column(6, numericInput("dose", "Dose (mg)", 300, min = 0, step = 50)),
                 column(6, numericInput("n", "Doses", 1, min = 1, max = 26, step = 1))),
        numericInput("every", "Every (days)", 28, min = 1, max = 84, step = 1, width = "100%"),
        muted("The trial: 30-300 mg SC and 300-750 mg IV single doses; 150-600 mg SC every 4 weeks x 3.")
      )
    } else {
      default <- switch(input$case, tce = 0.1, affinity = 0.1, cis = 0.1)
      tagList(
        fluidRow(column(6, numericInput("dose", "Dose (mg/kg IV)", default, min = 0, step = 0.05)),
                 column(6, numericInput("n", "Doses", if (input$case == "tce") 3 else 1, min = 1, max = 26, step = 1))),
        numericInput("every", "Every (days)", 7, min = 1, max = 84, step = 1, width = "100%")
      )
    }
  })

  # a case study's published parameters fill the setup tab
  observeEvent(input$case, {
    P <- case_params(M$model, cs()$sheet)
    for (k in c("CLp", "MWab", "Kd_R1", "Kd_R2", "BW")) updateNumericInput(session, k, value = P[[k]])
    updateNumericInput(session, "kon", value = P$kon_R1)
    updateNumericInput(session, "days", value = if (input$case == "soluble") 84 else 28)
  })

  params <- reactive({
    shiny::req(input$CLp, input$Kd_R1, input$Kd_R2, input$kon, input$BW > 0)
    kd3 <- if (input$case == "affinity" && !is.null(input$cd3)) {
      if (input$cd3 == "none") NA else as.numeric(input$cd3)
    } else 40
    P <- case_setup(M$model, input$case, cd3_kd = kd3,
                    chi_e = if (is.null(input$chi)) 1000 else as.numeric(input$chi),
                    format = if (is.null(input$format)) "bispecific" else input$format)
    P$CLp <- input$CLp; P$MWab <- input$MWab; P$BW <- input$BW
    # the sidebar's molecule choices (CD3 affinity, bivalent format) set the
    # affinities they are about; otherwise the setup tab does
    if (!(input$case == "cis" && identical(input$format, "bivalent"))) {
      P$Kd_R1 <- input$Kd_R1
      if (input$case != "affinity") P$Kd_R2 <- input$Kd_R2
    }
    P$kon_R1 <- input$kon
    if (P$kon_R2 > 0) P$kon_R2 <- input$kon
    P
  }) |> debounce(600)

  regimen <- reactive({
    shiny::req(input$dose >= 0, input$n >= 1, input$every >= 1)
    list(unit = if (input$case == "soluble") "mg" else "mg/kg",
         route = if (input$case == "soluble" && !is.null(input$route)) input$route else "IV",
         dose = input$dose, n = input$n, every = input$every, days = input$days)
  }) |> debounce(600)

  sim <- reactive({
    P <- params(); rg <- regimen()
    ev <- dose_events(P, rg$dose, rg$unit, rg$route, rg$n, rg$every)
    times <- c(0.01, seq(0.25, rg$days, by = 0.25))
    r <- withProgress(message = "Simulating", value = 0.4,
                      M$engine$solve(M$model, P, ev, times, rtol = case_tol(input$case)$rtol, atol = case_tol(input$case)$atol, nonneg = TRUE))
    list(r = r, t = times, P = P, rg = rg)
  })

  output$engine_note <- renderText(sprintf("  (%s)", M$engine$name))
  output$case_note <- renderUI(tags$p(tags$strong(cs()$species, ". "), cs()$note))

  st <- function(x, nm) x$r$states[, nm, 1]
  last <- function(v) v[length(v)]

  blood_drug <- function(x) if (cs()$kind == "soluble") x$r$total_D1_cen_ugml[, 1] else x$r$D1_cen_ugml[, 1]

  output$vb1 <- renderValueBox({
    x <- sim(); stat_box(signif(max(blood_drug(x)), 3), " ug/mL", "Peak in blood",
                         if (cs()$kind == "soluble") "Total drug" else "Free drug", "arrow-up")
  })
  output$vb2 <- renderValueBox({
    x <- sim(); stat_box(signif(last(blood_drug(x)), 3), " ug/mL", "Blood at the end", sprintf("Day %g", x$rg$days), "arrow-down")
  })
  output$vb3 <- renderValueBox({
    x <- sim(); k <- cs()$kind
    if (k == "soluble") stat_box(sprintf("%.0f", 100 - last(x$r$TN_sR2_cen[, 1])), "%", "IL-17AA neutralised", "Blood, at the end", "bullseye")
    else if (k == "trans") stat_box(sprintf("%.0f", last(x$r$RO_mR2_cen[, 1])), "%", "CD3 occupied", "Blood, at the end", "bullseye")
    else stat_box(sprintf("%.0f", x$r$RO_mR1_eff[which.min(abs(x$t - 21)), 1]), "%", "Therapeutic receptor, tumour", "Occupancy on day 21", "bullseye")
  })
  output$vb4 <- renderValueBox({
    x <- sim(); k <- cs()$kind
    if (k == "soluble") stat_box(sprintf("%.0f", 100 - last(x$r$TN_sR1_cen[, 1])), "%", "IL-13 neutralised", "Blood, at the end", "bullseye")
    else if (k == "trans") {
      tot <- last(x$r$CL_D1_tot_nmol[, 1])
      stat_box(sprintf("%.0f", 100 * last(st(x, "CL_D1_TMDD_R2_cen_nmol")) / max(tot, 1e-30)), "%", "Cleared through CD3",
               "Share of all drug cleared so far", "percent")
    } else stat_box(sprintf("%.0f", x$r$RO_mR1_cen[which.min(abs(x$t - 21)), 1]), "%", "Therapeutic receptor, blood", "Occupancy on day 21", "bullseye")
  })

  output$pk_plot <- renderPlot({
    x <- sim(); k <- cs()$kind
    lv <- c(if (k == "soluble") "Blood (total)" else "Blood (free)", "Leaky tissue", "Tight tissue", if (k == "cis") "Tumour")
    d <- rbind(data.frame(time = x$t, value = blood_drug(x), series = lv[1]),
               data.frame(time = x$t, value = x$r$D1_lea_ugml[, 1], series = lv[2]),
               data.frame(time = x$t, value = x$r$D1_tig_ugml[, 1], series = lv[3]),
               if (k == "cis") data.frame(time = x$t, value = x$r$D1_eff_ugml[, 1], series = lv[4]))
    d$value <- ifelse(d$value > max(d$value) * 1e-6, d$value, NA)
    series_plot(d, "Day", "Drug (ug/mL, log scale)", lv) + scale_y_log10(labels = plain_number)
  })
  output$pk_note <- renderUI({
    muted(switch(cs()$kind,
      soluble = "Total drug in blood (free plus bound to cytokines) is what the trial measured; the tissue lines are free drug.",
      trans = "Free drug. CD3 and CD20 were placed in blood only, as in the paper's calibration.",
      cis = "Free drug. The tumour is the model's efficacy space (4 mL/kg, 7% of lymph flow)."))
  })

  output$te_plot <- renderPlot({
    x <- sim(); k <- cs()$kind
    if (k == "soluble") {
      lv <- c("IL-13, blood", "IL-17AA, blood", "IL-13, leaky tissue", "IL-17AA, leaky tissue")
      d <- rbind(data.frame(time = x$t, value = x$r$TN_sR1_cen[, 1], series = lv[1]),
                 data.frame(time = x$t, value = x$r$TN_sR2_cen[, 1], series = lv[2]),
                 data.frame(time = x$t, value = x$r$TN_sR1_lea[, 1], series = lv[3]),
                 data.frame(time = x$t, value = x$r$TN_sR2_lea[, 1], series = lv[4]))
      series_plot(d, "Day", "Free cytokine (% of baseline)", lv) + scale_y_continuous(limits = c(0, NA))
    } else if (k == "trans") {
      lv <- c("CD20 occupied", "CD3 occupied")
      d <- rbind(data.frame(time = x$t, value = x$r$RO_mR1_cen[, 1], series = lv[1]),
                 data.frame(time = x$t, value = x$r$RO_mR2_cen[, 1], series = lv[2]))
      series_plot(d, "Day", "Receptors bound, blood (%)", lv) + scale_y_continuous(limits = c(0, 100))
    } else {
      lv <- c("Therapeutic receptor, tumour", "Therapeutic receptor, blood", "Targeting receptor, tumour")
      d <- rbind(data.frame(time = x$t, value = x$r$RO_mR1_eff[, 1], series = lv[1]),
                 data.frame(time = x$t, value = x$r$RO_mR1_cen[, 1], series = lv[2]),
                 data.frame(time = x$t, value = x$r$RO_mR2_eff[, 1], series = lv[3]))
      series_plot(d, "Day", "Receptor occupancy (%)", lv) + scale_y_continuous(limits = c(0, 100))
    }
  })
  output$te_note <- renderUI({
    muted(switch(cs()$kind,
      soluble = "Neutralisation = 100 - free/baseline. Total cytokine rises many-fold, held in drug complexes that clear like antibody.",
      trans = "Occupancy = fall of free receptor below its baseline. Drug bridging CD20 and CD3 is the trimer that makes a synapse.",
      cis = "The aim: occupy the therapeutic receptor in the tumour (above 60%) while sparing it in blood (below 40%), the paper's window."))
  })

  # --- dose response -------------------------------------------------------------

  output$dose_title <- renderText(switch(cs()$kind,
    soluble = "Target neutralisation on day 28 after single IV doses",
    trans = "Where the drug is cleared: share of clearance by day 21 after single IV doses",
    cis = "Therapeutic receptor occupancy on day 21 after single IV doses"))

  output$dose_plot <- renderPlot({
    P0 <- params()
    k <- cs()$kind
    if (k == "soluble") {
      doses <- c(30, 100, 300, 750, 1500)
      Pn <- P0[rep(1, length(doses)), ]
      ev <- do.call(rbind, lapply(seq_along(doses), function(i) dose_events(P0, doses[i], "mg", "IV", 1, 28, ID = i)))
      r <- withProgress(message = "Simulating the dose range", value = 0.4,
                        M$engine$solve(M$model, Pn, ev, c(0.01, 28), rtol = case_tol(input$case)$rtol, atol = case_tol(input$case)$atol, nonneg = TRUE))
      d <- rbind(data.frame(dose = doses, value = 100 - r$TN_sR1_cen[2, ], series = "IL-13, blood"),
                 data.frame(dose = doses, value = 100 - r$TN_sR2_cen[2, ], series = "IL-17AA, blood"),
                 data.frame(dose = doses, value = 100 - r$TN_sR1_lea[2, ], series = "IL-13, leaky tissue"),
                 data.frame(dose = doses, value = 100 - r$TN_sR2_lea[2, ], series = "IL-17AA, leaky tissue"),
                 data.frame(dose = doses, value = 100 - r$TN_sR2_tig[2, ], series = "IL-17AA, tight tissue"))
      d$time <- d$dose
      series_plot(d, "IV dose (mg, log scale)", "Neutralised on day 28 (%)", unique(d$series)) +
        geom_point(size = 2.2, na.rm = TRUE) + scale_x_log10(labels = plain_number) + scale_y_continuous(limits = c(0, 100))
    } else if (k == "trans") {
      doses <- c(0.01, 0.1, 1, 10)
      Pn <- P0[rep(1, length(doses)), ]
      ev <- do.call(rbind, lapply(seq_along(doses), function(i) dose_events(P0, doses[i], "mg/kg", "IV", 1, 7, ID = i)))
      r <- withProgress(message = "Simulating the dose range", value = 0.4,
                        M$engine$solve(M$model, Pn, ev, c(0.01, 21), rtol = case_tol(input$case)$rtol, atol = case_tol(input$case)$atol, nonneg = TRUE))
      s <- r$states[2, , ]
      tot <- s["CL_D1_NS_nmol", ] + s["CL_D1_TMDD_R1_cen_nmol", ] + s["CL_D1_TMDD_R2_cen_nmol", ]
      d <- rbind(data.frame(dose = doses, share = s["CL_D1_NS_nmol", ] / tot, path = "Nonspecific"),
                 data.frame(dose = doses, share = s["CL_D1_TMDD_R1_cen_nmol", ] / tot, path = "Through CD20 (B cells)"),
                 data.frame(dose = doses, share = s["CL_D1_TMDD_R2_cen_nmol", ] / tot, path = "Through CD3 (T cells)"))
      d$path <- factor(d$path, levels = c("Through CD3 (T cells)", "Through CD20 (B cells)", "Nonspecific"))
      d$dose_f <- factor(paste(plain_number(d$dose), "mg/kg"), levels = paste(plain_number(doses), "mg/kg"))
      ggplot(d, aes(dose_f, share, fill = path)) + geom_col(width = 0.6) +
        scale_fill_manual(values = c("Through CD3 (T cells)" = SERIES[2], "Through CD20 (B cells)" = SERIES[4], "Nonspecific" = PAL$ink_3)) +
        scale_y_continuous(labels = percent) + labs(x = "Single IV dose", y = "Share of drug cleared by day 21") + theme_sim(12)
    } else {
      doses <- c(0.01, 0.03, 0.1, 0.3, 1, 3, 10)
      runs <- lapply(c("bispecific", "bivalent"), function(fm) {
        P <- case_setup(M$model, "cis", chi_e = if (is.null(input$chi)) 1000 else as.numeric(input$chi), format = fm)
        P$CLp <- P0$CLp
        Pn <- P[rep(1, length(doses)), ]
        ev <- do.call(rbind, lapply(seq_along(doses), function(i) dose_events(P, doses[i], "mg/kg", "IV", 1, 7, ID = i)))
        r <- M$engine$solve(M$model, Pn, ev, c(0.01, 21), rtol = case_tol(input$case)$rtol, atol = case_tol(input$case)$atol, nonneg = TRUE)
        lab <- if (fm == "bispecific") "Bispecific" else "Bivalent"
        rbind(data.frame(time = doses, value = r$RO_mR1_eff[2, ], series = paste(lab, "- tumour")),
              data.frame(time = doses, value = r$RO_mR1_cen[2, ], series = paste(lab, "- blood")))
      })
      d <- do.call(rbind, runs)
      series_plot(d, "Single IV dose (mg/kg, log scale)", "Therapeutic receptor occupied, day 21 (%)", unique(d$series)) +
        annotate("rect", xmin = min(doses) / 1.5, xmax = max(doses) * 1.5, ymin = 60, ymax = 100, alpha = 0.06, fill = PAL$green) +
        geom_point(size = 2.2) + scale_x_log10(labels = plain_number) + scale_y_continuous(limits = c(0, 100))
    }
  })
  output$dose_note <- renderUI({
    muted(switch(cs()$kind,
      soluble = "The paper (Fig. 3d-e): 750 mg IV neutralises about 80% of IL-17AA and about 90% of IL-13 in serum on day 28. IL-17AA turns over in minutes, so it is harder to hold down than IL-13.",
      trans = "The paper (Fig. 6): with a moderate CD3 affinity most of the drug is cleared through T cells below 10 mg/kg; weaken CD3 binding (400 nM) and nonspecific clearance dominates. The Molecule/CD3 choice in the sidebar applies here.",
      cis = "The paper (Fig. 7): avidity barely changes occupancy in blood but controls it in the tumour; at high doses the tumour advantage is lost. A bivalent 10 nM antibody engages blood as much as tumour. Shaded: above 60% (the efficacy threshold)."))
  })

  output$par_table <- renderDT({
    P <- params()
    small_table(data.frame(Parameter = names(P), Value = signif(unlist(P[1, ]), 4)), page = 200)
  })
}

shinyApp(ui, server)
