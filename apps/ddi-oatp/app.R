# ============================================================================
# Rifampicin -> glibenclamide: OATP inhibition versus CYP2C9/3A induction
# (Asaumi et al. 2018, CPT Pharmacometrics Syst Pharmacol)
#
# Glibenclamide is taken up into hepatocytes by OATP transporters and
# metabolised by CYP2C9 and CYP3A. Rifampicin inhibits OATP while it is
# present (exposure up) and induces CYP2C9 and CYP3A over days (exposure
# down), so the net interaction depends entirely on timing - a classic
# "complex DDI" that static models cannot capture.
# ============================================================================

library(shiny)
library(shinydashboard)
library(DT)
library(ggplot2)
library(scales)

shared <- function(f) source(if (file.exists(file.path("R", f))) file.path("R", f) else file.path("..", "..", "shared", f))
for (f in c("theme.R", "ode_engine.R", "app_helpers.R", "rifampicin.R")) shared(f)

M <- load_model("ddi_rifampicin_glibenclamide.cpp")

# Asaumi 2018 Figure 6: glibenclamide 1.25 mg PO in each case.
SCENARIOS <- list(
  iv = list(label = "Single IV rifampicin with glibenclamide (OATP inhibition)", po_days = 0, iv = TRUE, gap_h = 0,
            paper = 2.08, observed = "2.18 +/- 1.09"),
  poiv = list(label = "6 days oral rifampicin, then IV with glibenclamide on day 7", po_days = 6, iv = TRUE, gap_h = 0,
              paper = 0.90, observed = "0.72 +/- 0.32"),
  po = list(label = "7 days oral rifampicin, glibenclamide on day 9 (induction)", po_days = 7, iv = FALSE, gap_h = 48,
            paper = 0.50, observed = "0.35 +/- 0.19")
)

ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = tags$div(tags$span("PBPK", style = "font-weight: bold; font-size: 24px;"),
                                   tags$span(" DDI: OATP + induction", style = "font-size: 18px;")),
                  titleWidth = 350),
  dashboardSidebar(
    width = 320, custom_css,
    sidebarMenu(
      menuItem("Interaction", tabName = "main", icon = icon("right-left")),
      menuItem("Timing matters", tabName = "timing", icon = icon("clock")),
      menuItem("About", tabName = "about", icon = icon("info-circle"))
    ),
    tags$div(
      style = "padding: 18px; padding-top: 8px;",
      tags$h4("Scenario"),
      radioButtons("scen", NULL, choices = c(stats::setNames(names(SCENARIOS), vapply(SCENARIOS, `[[`, "", "label")),
                                             "Custom" = "custom"), selected = "iv"),
      conditionalPanel(
        "input.scen == 'custom'",
        numericInput("rif_dose", "Rifampicin dose (mg)", 600, min = 0, max = 1200, step = 50),
        sliderInput("po_days", "Days of oral rifampicin first", min = 0, max = 14, value = 7, step = 1, width = "100%"),
        checkboxInput("iv_with", "Plus an IV rifampicin dose given with glibenclamide", FALSE),
        sliderInput("gap_h", "Hours from the last oral dose to glibenclamide", min = 0, max = 168, value = 48, step = 6, width = "100%")
      ),
      tags$h4("Victim: glibenclamide"),
      numericInput("glb_dose", "Oral dose (mg)", 1.25, min = 0.25, max = 10, step = 0.25, width = "100%"),
      numericInput("bw", "Body weight (kg)", 70, min = 40, max = 150, step = 5, width = "100%")
    )
  ),
  dashboardBody(
    custom_css,
    tabItems(
      tabItem(
        tabName = "main",
        fluidRow(valueBoxOutput("aucr", width = 3), valueBoxOutput("paper", width = 3),
                 valueBoxOutput("oatp", width = 3), valueBoxOutput("cyp", width = 3)),
        fluidRow(
          box(title = tags$div(tags$strong("Glibenclamide in blood"),
                               tags$span(textOutput("engine_note", inline = TRUE),
                                         style = "color: #6b7078; font-size: 13px; font-weight: normal;")),
              status = "primary", solidHeader = TRUE, width = 7, plotOutput("glb_plot", height = "380px")),
          box(title = "What rifampicin is doing", status = "primary", solidHeader = TRUE, width = 5,
              plotOutput("mech_plot", height = "380px"),
              muted("OATP activity = 1 / (1 + unbound rifampicin / Ki,u) at the liver inlet; enzyme activities from the turnover model."))
        )
      ),
      tabItem(
        tabName = "timing",
        box(title = "Glibenclamide AUC ratio by when it is taken", status = "primary", solidHeader = TRUE, width = 12,
            actionButton("run_t", "Run timing scan", icon = icon("play")),
            muted("Seven days of oral rifampicin (the dose set under Custom, or 600 mg), then glibenclamide at different times after the last rifampicin dose. Taken with rifampicin still present, OATP inhibition partly offsets the induction; taken a day or more later, induction dominates; after a week or two, both have worn off."),
            plotOutput("t_plot", height = "400px"))
      ),
      about_tab(
        "Rifampicin -> glibenclamide (OATP inhibition and CYP induction)",
        "A complex drug-drug interaction in which one perpetrator has two opposite effects
         with different time courses. The PBPK model resolves them mechanistically:
         competitive inhibition of OATP-mediated hepatic uptake by unbound rifampicin in
         each liver zone, and induction of CYP2C9 and CYP3A through enzyme turnover.",
        list("Rifampicin: as in the CYP3A-induction app, plus CYP2C9 induction (Emax 2.41, kdeg 0.00666 /h)",
             "Glibenclamide: blood, 5 hepatic extracellular + 5 hepatocyte compartments, muscle, skin, adipose, portal vein, gut lumen",
             "Hepatic uptake: OATP (PSact,inf) competitively inhibited by unbound rifampicin (Ki,u 0.226 uM) plus passive diffusion; efflux; metabolism fm CYP2C9 0.85 / CYP3A 0.15 (beta 0.2)",
             "Gut: Qgut model with CYP3A-dependent intestinal metabolism"),
        tags$span("Asaumi R, Toshimoto K, Tobe Y, et al. Comprehensive PBPK model of rifampicin for
                   quantitative prediction of complex drug-drug interactions: CYP3A/2C9 induction and
                   OATP inhibition effects. CPT Pharmacometrics Syst Pharmacol 2018;7:186-196.
                   Glibenclamide transcribed from the Supplementary Text model code; OATP Ki,u from
                   Yoshikado et al. 2016 as quoted in the paper."),
        M$engine$name
      )
    )
  )
)

server <- function(input, output, session) {

  run <- function(rif_dose, po_days, iv_with, gap_h, glb_dose, bw) {
    t_glb <- if (po_days > 0) (po_days - 1) * 24 + gap_h else 0
    if (po_days > 0 && iv_with) t_glb <- po_days * 24
    ev <- rbind(rif_events(rif_dose, po_days, "po", 0, bw),
                if (iv_with) rif_events(rif_dose, 1, "iv", t_glb, bw),
                data.frame(time = t_glb + GLB_TLAG, cmt = "GLB_GUT", amt = glb_dose, rate = 0))
    tt <- ddi_times(t_glb, 48)
    r <- M$engine$solve(M$model, data.frame(BW = bw), ev, tt, rtol = 1e-5, atol = 1e-10)
    w <- tt >= t_glb
    list(r = r, t = tt, w = w, t_glb = t_glb, auc = auc_trap(tt[w] - t_glb, r$GLB[w, 1]))
  }

  setup <- reactive({
    shiny::req(input$glb_dose > 0, input$bw > 0)
    if (input$scen == "custom") {
      list(rif = input$rif_dose, po = input$po_days, iv = isTRUE(input$iv_with), gap = input$gap_h, paper = NA, obs = NA)
    } else {
      s <- SCENARIOS[[input$scen]]
      list(rif = 600, po = s$po_days, iv = s$iv, gap = s$gap_h, paper = s$paper, obs = s$observed)
    }
  }) |> debounce(500)

  sim <- reactive({
    s <- setup()
    withProgress(message = "Simulating rifampicin and glibenclamide", value = 0.3, {
      ddi <- run(s$rif, s$po, s$iv, s$gap, input$glb_dose, input$bw)
      ctl <- run(0, 0, FALSE, 0, input$glb_dose, input$bw)
    })
    list(ddi = ddi, ctl = ctl, s = s)
  })

  output$engine_note <- renderText(sprintf("  (%s)", M$engine$name))

  output$aucr <- renderValueBox({
    x <- sim(); v <- x$ddi$auc / x$ctl$auc
    stat_box(sprintf("%.2f", v), "", "Glibenclamide AUC ratio", "With vs without rifampicin (0-48 h)", "chart-area",
             color = if (v > 1.25) "red" else if (v < 0.8) "yellow" else "blue")
  })
  output$paper <- renderValueBox({
    x <- sim()
    if (is.na(x$s$paper)) return(stat_box("-", "", "Published value", "Custom scenario", "book"))
    stat_box(sprintf("%.2f", x$s$paper), "", "Paper's simulated AUCR", paste("Observed:", x$s$obs), "book")
  })
  output$oatp <- renderValueBox({
    x <- sim()
    k <- x$ddi$w & x$ddi$t <= x$ddi$t_glb + 12
    stat_box(sprintf("%.0f", 100 * min(x$ddi$r$OATP_ACTIVITY[k, 1])), "%", "Lowest OATP activity",
             "In the 12 h after glibenclamide", "door-open")
  })
  output$cyp <- renderValueBox({
    x <- sim(); k <- which(x$ddi$w)[1]
    stat_box(sprintf("%.2f", x$ddi$r$CYP2C9_LIVER[k, 1]), "x", "Hepatic CYP2C9", "At the glibenclamide dose", "vial")
  })

  output$glb_plot <- renderPlot({
    x <- sim()
    lv <- c("With rifampicin", "Without")
    d <- rbind(data.frame(time = x$ddi$t[x$ddi$w] - x$ddi$t_glb, value = x$ddi$r$GLB[x$ddi$w, 1], series = lv[1]),
               data.frame(time = x$ctl$t[x$ctl$w] - x$ctl$t_glb, value = x$ctl$r$GLB[x$ctl$w, 1], series = lv[2]))
    series_plot(d, "Hours after the glibenclamide dose", "Glibenclamide (ng/mL)", lv) +
      scale_x_continuous(breaks = seq(0, 48, 8))
  })

  output$mech_plot <- renderPlot({
    x <- sim(); r <- x$ddi$r
    lv <- c("OATP activity", "CYP2C9 (liver)", "CYP3A (liver)")
    d <- rbind(data.frame(time = x$ddi$t / 24, value = r$OATP_ACTIVITY[, 1], series = lv[1]),
               data.frame(time = x$ddi$t / 24, value = r$CYP2C9_LIVER[, 1], series = lv[2]),
               data.frame(time = x$ddi$t / 24, value = r$CYP3A_LIVER[, 1], series = lv[3]))
    series_plot(d, "Day", "Fraction of baseline", lv) +
      geom_hline(yintercept = 1, colour = PAL$ink_3) +
      geom_vline(xintercept = x$ddi$t_glb / 24, colour = PAL$ink_3, linetype = "22") +
      labs(subtitle = "Dashed line: glibenclamide dose")
  })

  timing <- eventReactive(input$run_t, {
    dose <- if (input$scen == "custom") input$rif_dose else 600
    gaps <- c(0, 6, 12, 24, 48, 96, 168, 336)
    ctl <- run(0, 0, FALSE, 0, input$glb_dose, input$bw)
    withProgress(message = "Scanning glibenclamide timing", value = 0, {
      do.call(rbind, lapply(gaps, function(g) {
        incProgress(1 / length(gaps), detail = sprintf("%d h", g))
        x <- run(dose, 7, FALSE, g, input$glb_dose, input$bw)
        data.frame(gap = g, aucr = x$auc / ctl$auc)
      }))
    })
  })

  output$t_plot <- renderPlot({
    d <- timing()
    ggplot(d, aes(gap / 24, aucr)) +
      geom_hline(yintercept = 1, colour = PAL$ink_3, linetype = "22") +
      geom_line(colour = PAL$blue_ink, linewidth = 1.1) + geom_point(colour = PAL$blue_ink, size = 3) +
      labs(x = "Days from the last rifampicin dose to glibenclamide", y = "Glibenclamide AUC ratio") +
      theme_sim(12)
  })
}

shinyApp(ui, server)
