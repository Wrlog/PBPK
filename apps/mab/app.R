# ============================================================================
# Monoclonal antibodies: second-generation minimal PBPK with TMDD
# (Cao, Balthasar & Jusko 2013; Cao & Jusko 2014)
#
# Antibody leaves plasma by convection through vascular pores into the
# interstitial fluid of "tight" and "leaky" tissues and returns through
# lymph. Nonspecific clearance is from plasma. With a target in the
# interstitial fluid, binding and internalisation add a saturable clearance
# - the source of the nonlinear PK seen at low doses of many antibodies.
# ============================================================================

library(shiny)
library(shinydashboard)
library(DT)
library(ggplot2)
library(scales)

shared <- function(f) source(if (file.exists(file.path("R", f))) file.path("R", f) else file.path("..", "..", "shared", f))
for (f in c("theme.R", "ode_engine.R", "app_helpers.R")) shared(f)

M <- load_model("mab_mpbpk_tmdd.cpp")

# Cao & Jusko 2014, Table 2 (pTMDD, patients); Cao et al. 2013, Table 2
# (Model A, linear, 70 kg; CLp converted from L/h to mL/h/kg).
MABS <- list(
  "Trastuzumab (HER2, with TMDD)" = list(SIGMA1 = 0.95, SIGMA2 = 0.512, CLp_kg = 0.092, TMDD = 1,
                                          KSS = 0.99, KSYN = 0.376, KDEG = 0.0117, KINT = 0.0117),
  "Adecatumumab" = list(SIGMA1 = 0.883, SIGMA2 = 0.524, CLp_kg = 0.0300 * 1000 / 70),
  "Mepolizumab" = list(SIGMA1 = 0.950, SIGMA2 = 0.750, CLp_kg = 0.00851 * 1000 / 70),
  "Gevokizumab" = list(SIGMA1 = 0.931, SIGMA2 = 0.837, CLp_kg = 0.00668 * 1000 / 70),
  "GNbAC1" = list(SIGMA1 = 0.915, SIGMA2 = 0.831, CLp_kg = 0.00714 * 1000 / 70),
  "MEDI-528" = list(SIGMA1 = 0.987, SIGMA2 = 0.754, CLp_kg = 0.00537 * 1000 / 70),
  "Tefibazumab" = list(SIGMA1 = 0.902, SIGMA2 = 0.815, CLp_kg = 0.00933 * 1000 / 70),
  "PAmAb" = list(SIGMA1 = 0.950, SIGMA2 = 0.779, CLp_kg = 0.00867 * 1000 / 70),
  "PRO95780" = list(SIGMA1 = 0.984, SIGMA2 = 0.638, CLp_kg = 0.0124 * 1000 / 70),
  "Siltuximab" = list(SIGMA1 = 0.964, SIGMA2 = 0.673, CLp_kg = 0.0115 * 1000 / 70),
  "Visilizumab" = list(SIGMA1 = 0.949, SIGMA2 = 0.834, CLp_kg = 0.0152 * 1000 / 70)
)

ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = tags$div(tags$span("PBPK", style = "font-weight: bold; font-size: 24px;"),
                                   tags$span(" Antibody mPBPK", style = "font-size: 18px;")),
                  titleWidth = 350),
  dashboardSidebar(
    width = 320, custom_css,
    sidebarMenu(
      menuItem("Plasma and tissue", tabName = "main", icon = icon("vial")),
      menuItem("Nonlinearity", tabName = "nonlin", icon = icon("chart-line")),
      menuItem("Model setup", tabName = "setup", icon = icon("sliders")),
      menuItem("About", tabName = "about", icon = icon("info-circle"))
    ),
    tags$div(
      style = "padding: 18px; padding-top: 8px;",
      tags$h4("Antibody"),
      selectInput("mab", NULL, choices = names(MABS), selected = names(MABS)[1], width = "100%"),
      checkboxInput("tmdd", "Target-mediated disposition in the ISF", TRUE),
      tags$h4("Regimen (IV infusion)"),
      fluidRow(
        column(6, numericInput("load", "Loading dose (mg/kg)", 4, min = 0, max = 30, step = 0.5)),
        column(6, numericInput("maint", "Maintenance (mg/kg)", 2, min = 0, max = 30, step = 0.5))
      ),
      fluidRow(
        column(6, numericInput("ii", "Every (weeks)", 1, min = 1, max = 8, step = 1)),
        column(6, numericInput("n", "Doses", 12, min = 1, max = 52, step = 1))
      ),
      numericInput("bw", "Body weight (kg)", 70, min = 30, max = 150, step = 5, width = "100%"),
      muted("The loading dose is the first infusion; the maintenance dose follows every interval.")
    )
  ),
  dashboardBody(
    custom_css,
    tabItems(
      tabItem(
        tabName = "main",
        fluidRow(valueBoxOutput("cmax", width = 3), valueBoxOutput("trough", width = 3),
                 valueBoxOutput("isf", width = 3), valueBoxOutput("ro", width = 3)),
        fluidRow(
          box(title = tags$div(tags$strong("Antibody in plasma and interstitial fluid"),
                               tags$span(textOutput("engine_note", inline = TRUE),
                                         style = "color: #6b7078; font-size: 13px; font-weight: normal;")),
              status = "primary", solidHeader = TRUE, width = 7, plotOutput("conc_plot", height = "400px"),
              muted("Free antibody. Leaky tissues (liver, kidney, heart, ...) equilibrate much faster and higher than tight ones (muscle, skin, fat, brain), whose vessels reflect most of the antibody.")),
          box(title = "Target engagement in the interstitial fluid", status = "primary", solidHeader = TRUE, width = 5,
              plotOutput("ro_plot", height = "400px"))
        )
      ),
      tabItem(
        tabName = "nonlin",
        box(title = "Dose-normalised plasma profiles after a single dose", status = "primary", solidHeader = TRUE, width = 12,
            plotOutput("nl_plot", height = "440px"),
            muted("Linear kinetics would make these curves identical. Target-mediated clearance is saturable, so low doses are cleared faster: the curves separate once the concentration falls to the level of the target."))
      ),
      tabItem(
        tabName = "setup",
        fluidRow(
          box(title = "Distribution and clearance", status = "primary", solidHeader = TRUE, width = 6,
              fluidRow(column(6, numericInput("SIGMA1", "Reflection, tight (sigma1)", 0.95, min = 0, max = 1, step = 0.01)),
                       column(6, numericInput("SIGMA2", "Reflection, leaky (sigma2)", 0.512, min = 0, max = 1, step = 0.01))),
              fluidRow(column(6, numericInput("CLp_kg", "Plasma clearance (mL/h/kg)", 0.092, min = 0, step = 0.01)),
                       column(6, numericInput("Kp", "Available ISF fraction (Kp)", 0.8, min = 0.1, max = 1, step = 0.1))),
              muted("Physiology for 70 kg: plasma 2.6 L, ISF 15.6 L (65% tight / 35% leaky), lymph 5.2 L, lymph flow 2.9 L/day (33% / 67%); lymph reflection 0.2. All scale with body weight.")),
          box(title = "Target (used when TMDD is on)", status = "primary", solidHeader = TRUE, width = 6,
              fluidRow(column(6, numericInput("KSS", "Kss (nM)", 0.99, min = 0.001, step = 0.1)),
                       column(6, numericInput("KSYN", "ksyn (nM/h)", 0.376, min = 0, step = 0.05))),
              fluidRow(column(6, numericInput("KDEG", "kdeg (1/h)", 0.0117, min = 0.0001, step = 0.001)),
                       column(6, numericInput("KINT", "kint (1/h)", 0.0117, min = 0, step = 0.001))),
              muted("Baseline target R0 = ksyn / kdeg in each ISF space. Quasi-steady-state binding (Gibiansky)."))
        )
      ),
      about_tab(
        "Second-generation minimal PBPK for monoclonal antibodies",
        "A minimal PBPK model that keeps the physiology that matters for IgG antibodies -
         plasma, lymph and the interstitial fluid of two lumped tissue groups, joined by
         convective flow through vascular pores - and fits plasma data with only three
         parameters. The 2014 extension adds target-mediated disposition in the
         interstitial fluid, where most membrane targets are.",
        list("Plasma: nonspecific clearance CLp; convective loss to ISF (1 - sigma) x lymph flow",
             "ISF: tight tissues (continuous endothelium; sigma1) and leaky tissues (fenestrated/discontinuous; sigma2), 65% / 35% of available ISF",
             "Lymph: collects ISF (reflection 0.2) and returns it to plasma",
             "TMDD (optional): target synthesis ksyn, degradation kdeg, antibody binding (Kss, quasi-steady state), complex internalisation kint, in both ISF spaces; only free antibody drains to lymph",
             "Presets: trastuzumab with its HER2 TMDD parameters (Cao & Jusko 2014), and ten linear antibodies (Cao et al. 2013, Model A)"),
        tags$span("Cao Y, Balthasar JP, Jusko WJ. Second-generation minimal physiologically-based
                   pharmacokinetic model for monoclonal antibodies. J Pharmacokinet Pharmacodyn
                   2013;40:597-607. Cao Y, Jusko WJ. Incorporating target-mediated drug disposition in a
                   minimal physiologically-based pharmacokinetic model for monoclonal antibodies.
                   J Pharmacokinet Pharmacodyn 2014;41:375-387."),
        M$engine$name
      )
    )
  )
)

server <- function(input, output, session) {

  observeEvent(input$mab, {
    x <- MABS[[input$mab]]
    for (k in c("SIGMA1", "SIGMA2", "CLp_kg")) updateNumericInput(session, k, value = x[[k]])
    has_t <- !is.null(x$TMDD)
    updateCheckboxInput(session, "tmdd", value = has_t)
    if (has_t) for (k in c("KSS", "KSYN", "KDEG", "KINT")) updateNumericInput(session, k, value = x[[k]])
  })

  params <- reactive({
    shiny::req(input$SIGMA1, input$SIGMA2, input$CLp_kg, input$bw > 0)
    data.frame(BW = input$bw, SIGMA1 = input$SIGMA1, SIGMA2 = input$SIGMA2, CLp_kg = input$CLp_kg, Kp = input$Kp,
               TMDD = as.numeric(isTRUE(input$tmdd)), KSS = input$KSS, KSYN = input$KSYN,
               KDEG = input$KDEG, KINT = input$KINT)
  }) |> debounce(500)

  nmol <- function(mg) mg * 1e6 / M$model$param[["MW"]]

  regimen <- reactive({
    shiny::req(input$n >= 1, input$ii >= 1)
    bw <- input$bw
    ev <- rbind(
      if (input$load > 0) data.frame(time = 0, cmt = "CENT", amt = nmol(input$load * bw), rate = nmol(input$load * bw) / 1.5, ii = 0, addl = 0),
      if (input$maint > 0 && input$n > 1) data.frame(time = input$ii * 168, cmt = "CENT", amt = nmol(input$maint * bw),
                                                     rate = nmol(input$maint * bw) / 0.5, ii = input$ii * 168, addl = input$n - 2)
    )
    list(ev = ev, end = input$ii * 168 * input$n + 168 * 2, last = (input$n - 1) * input$ii * 168)
  }) |> debounce(500)

  sim <- reactive({
    rg <- regimen()
    times <- seq(0, rg$end, by = 6)
    r <- withProgress(message = "Simulating antibody disposition", value = 0.4,
                      M$engine$solve(M$model, params(), rg$ev, times, rtol = 1e-5, atol = 1e-9, nonneg = TRUE))
    list(r = r, t = times, rg = rg)
  })

  output$engine_note <- renderText(sprintf("  (%s)", M$engine$name))

  trough_idx <- function(x) max(which(x$t <= x$rg$last - 1e-6), 1)

  output$cmax <- renderValueBox({
    x <- sim(); stat_box(sprintf("%.0f", max(x$r$CP_UGML[, 1])), " ug/mL", "Peak plasma", "Whole course", "arrow-up")
  })
  output$trough <- renderValueBox({
    x <- sim(); stat_box(sprintf("%.1f", x$r$CP_UGML[trough_idx(x), 1]), " ug/mL", "Trough before last dose", "Plasma", "arrow-down")
  })
  output$isf <- renderValueBox({
    x <- sim(); k <- trough_idx(x)
    stat_box(sprintf("%.0f", 100 * x$r$ISF_LEAKY_UGML[k, 1] / x$r$CP_UGML[k, 1]), "%", "Leaky ISF / plasma", "Free antibody, at trough", "percent")
  })
  output$ro <- renderValueBox({
    x <- sim(); k <- trough_idx(x)
    if (!isTRUE(input$tmdd)) return(stat_box("-", "", "Target occupancy", "TMDD off", "bullseye"))
    stat_box(sprintf("%.0f", 100 * x$r$RO_TIGHT[k, 1]), "%", "Target occupancy, tight ISF", "At trough (lowest)", "bullseye")
  })

  output$conc_plot <- renderPlot({
    x <- sim()
    lv <- c("Plasma", "Leaky ISF (free)", "Tight ISF (free)")
    d <- rbind(data.frame(time = x$t / 168, value = x$r$CP_UGML[, 1], series = lv[1]),
               data.frame(time = x$t / 168, value = x$r$ISF_LEAKY_UGML[, 1], series = lv[2]),
               data.frame(time = x$t / 168, value = x$r$ISF_TIGHT_UGML[, 1], series = lv[3]))
    d$value <- ifelse(d$value > max(d$value) * 1e-5, d$value, NA)
    series_plot(d, "Week", "Antibody (ug/mL, log scale)", lv) + scale_y_log10(labels = plain_number)
  })

  output$ro_plot <- renderPlot({
    x <- sim()
    if (!isTRUE(input$tmdd)) {
      return(ggplot() + annotate("text", x = 0, y = 0, label = "Target-mediated disposition is off", colour = PAL$ink_3) + theme_void())
    }
    lv <- c("Leaky ISF", "Tight ISF")
    d <- rbind(data.frame(time = x$t / 168, value = 100 * x$r$RO_LEAKY[, 1], series = lv[1]),
               data.frame(time = x$t / 168, value = 100 * x$r$RO_TIGHT[, 1], series = lv[2]))
    series_plot(d, "Week", "Target occupancy (%)", lv) + scale_y_continuous(limits = c(0, 100))
  })

  output$nl_plot <- renderPlot({
    P <- params()
    doses <- c(0.1, 0.3, 1, 3, 10)
    times <- seq(0, 168 * 6, by = 6)
    Pn <- P[rep(1, length(doses)), ]
    ev <- data.frame(time = 0, cmt = "CENT", amt = nmol(doses * P$BW), rate = nmol(doses * P$BW) / 1.5,
                     ID = seq_along(doses))
    r <- M$engine$solve(M$model, Pn, ev, times, rtol = 1e-5, atol = 1e-10, nonneg = TRUE)
    lv <- paste(plain_number(doses), "mg/kg")
    d <- do.call(rbind, lapply(seq_along(doses), function(k) {
      data.frame(time = times / 168, value = r$CP_UGML[, k] / doses[k], series = lv[k])
    }))
    d$value <- ifelse(d$value > max(d$value) * 1e-5, d$value, NA)
    series_plot(d, "Week", "Plasma per mg/kg (ug/mL per mg/kg, log)", lv) + scale_y_log10(labels = plain_number)
  })
}

shinyApp(ui, server)
