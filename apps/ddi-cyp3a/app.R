# ============================================================================
# Rifampicin -> midazolam: CYP3A induction in liver and gut
# (Asaumi et al. 2018, CPT Pharmacometrics Syst Pharmacol)
#
# Rifampicin, the reference CYP3A inducer, raises CYP3A in hepatocytes and
# enterocytes through an enzyme-turnover model, and also induces its own
# metabolism (UGT). Midazolam, the reference CYP3A substrate, shows the
# consequence: after a week of rifampicin an oral dose gives a small
# fraction of its usual exposure, and the effect wears off over days as the
# extra enzyme is degraded.
# ============================================================================

library(shiny)
library(shinydashboard)
library(DT)
library(ggplot2)
library(scales)

shared <- function(f) source(if (file.exists(file.path("R", f))) file.path("R", f) else file.path("..", "..", "shared", f))
for (f in c("theme.R", "ode_engine.R", "app_helpers.R", "rifampicin.R")) shared(f)

M <- load_model("ddi_rifampicin_midazolam.cpp")

ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = tags$div(tags$span("PBPK", style = "font-weight: bold; font-size: 24px;"),
                                   tags$span(" DDI: CYP3A induction", style = "font-size: 18px;")),
                  titleWidth = 350),
  dashboardSidebar(
    width = 320, custom_css,
    sidebarMenu(
      menuItem("Interaction", tabName = "main", icon = icon("right-left")),
      menuItem("Dose-response", tabName = "dr", icon = icon("chart-line")),
      menuItem("About", tabName = "about", icon = icon("info-circle"))
    ),
    tags$div(
      style = "padding: 18px; padding-top: 8px;",
      tags$h4("Perpetrator: rifampicin"),
      sliderInput("rif_dose", "Oral dose, once daily (mg)", min = 0, max = 600, value = 600, step = 5, width = "100%"),
      sliderInput("rif_days", "Days of rifampicin", min = 1, max = 14, value = 7, step = 1, width = "100%"),
      tags$h4("Victim: midazolam"),
      radioButtons("mdz_route", NULL, choices = c("Oral" = "po", "IV bolus" = "iv"), selected = "po", inline = TRUE),
      numericInput("mdz_dose", "Dose (mg)", value = 3, min = 0.1, max = 15, step = 0.5, width = "100%"),
      sliderInput("washout", "Given this many days after the last rifampicin dose", min = 0, max = 14, value = 1, step = 1, width = "100%"),
      muted("Day 1 means the morning after the last rifampicin dose. Longer gaps show the induction wearing off."),
      numericInput("bw", "Body weight (kg)", value = 70, min = 40, max = 150, step = 5, width = "100%")
    )
  ),
  dashboardBody(
    custom_css,
    tabItems(
      tabItem(
        tabName = "main",
        fluidRow(valueBoxOutput("aucr", width = 3), valueBoxOutput("cmaxr", width = 3),
                 valueBoxOutput("liver", width = 3), valueBoxOutput("gut", width = 3)),
        fluidRow(
          box(title = tags$div(tags$strong("Midazolam in blood"),
                               tags$span(textOutput("engine_note", inline = TRUE),
                                         style = "color: #6b7078; font-size: 13px; font-weight: normal;")),
              status = "primary", solidHeader = TRUE, width = 7,
              checkboxInput("log", "Log concentration axis", TRUE),
              plotOutput("mdz_plot", height = "360px")),
          box(title = "CYP3A activity (fold of baseline)", status = "primary", solidHeader = TRUE, width = 5,
              plotOutput("cyp_plot", height = "400px"),
              muted("Enzyme-turnover model: synthesis is induced by unbound rifampicin; the extra enzyme is lost with the natural degradation half-life (about 44 h in liver, 24 h in gut)."))
        ),
        fluidRow(
          box(title = "Rifampicin in blood", status = "primary", solidHeader = TRUE, width = 7,
              plotOutput("rif_plot", height = "280px"),
              muted("Rifampicin also induces the UGT that clears it (auto-induction), so its own exposure falls over the first days.")),
          box(title = "Exposure", status = "info", solidHeader = TRUE, width = 5, DT::dataTableOutput("tbl"))
        )
      ),
      tabItem(
        tabName = "dr",
        box(title = "Midazolam AUC ratio across rifampicin doses", status = "primary", solidHeader = TRUE, width = 12,
            actionButton("run_dr", "Run dose-response", icon = icon("play")),
            muted("Seven rifampicin doses from 5 to 600 mg, each with the regimen, midazolam route and timing set in the sidebar. The paper reported simulated oral-midazolam AUC ratios close to the observed ones across 5-600 mg; at 600 mg the observed oral AUC falls by about 90-95%."),
            plotOutput("dr_plot", height = "400px"), DT::dataTableOutput("dr_tbl"))
      ),
      about_tab(
        "Rifampicin -> midazolam (CYP3A induction)",
        "A PBPK DDI model for the most-studied induction interaction. Rifampicin's own
         model has a dispersion-type liver, OATP-mediated uptake into hepatocytes, a
         segregated-flow gut and auto-induction of its UGT metabolism; its unbound
         concentrations in hepatocytes and enterocytes induce CYP3A through an
         enzyme-turnover model. Midazolam's intrinsic clearance in liver and gut scales
         with the local CYP3A level.",
        list("Rifampicin: blood, 5 hepatic extracellular + 5 hepatocyte compartments, muscle, skin, adipose, serosa, gut lumen, enterocytes and mucosal blood",
             "Hepatic uptake: saturable OATP transport plus passive diffusion (Km,u 146 ng/mL), efflux and UGT metabolism in hepatocytes",
             "Induction: dE/dt = kdeg (1 + Emax Cu / (Cu + EC50,u) - E); Emax 4.57 (CYP3A), 1.34 (UGT); EC50,u 52.6 ng/mL; kdeg 0.0158 /h liver, 0.0288 /h gut",
             "Midazolam: blood, 5 liver compartments, muscle, skin, adipose, portal vein, gut lumen (Qgut model); fm,CYP3A 0.93 in liver, 1 in gut",
             "All physiology per kg of body weight (Table S1)"),
        tags$span("Asaumi R, Toshimoto K, Tobe Y, et al. Comprehensive PBPK model of rifampicin for
                   quantitative prediction of complex drug-drug interactions: CYP3A/2C9 induction and
                   OATP inhibition effects. CPT Pharmacometrics Syst Pharmacol 2018;7:186-196.
                   Transcribed from the model code in the Supplementary Text, with parameters from
                   Supplementary Tables S1-S2 and Table 1."),
        M$engine$name
      )
    )
  )
)

server <- function(input, output, session) {

  run <- function(rif_dose, rif_days, route, mdz_dose, washout, bw) {
    t_last <- (rif_days - 1) * 24
    t_mdz <- if (rif_dose > 0) t_last + 24 * washout else 0
    victim <- if (route == "po") {
      data.frame(time = t_mdz + MDZ_TLAG, cmt = "MDZ_GUT", amt = mdz_dose, rate = 0)
    } else {
      data.frame(time = t_mdz, cmt = "MDZ_CENT", amt = mdz_dose / (MDZ_VC_PER_KG * bw), rate = 0)
    }
    ev <- rbind(rif_events(rif_dose, rif_days, "po", 0, bw), victim)
    tt <- ddi_times(t_mdz, 24)
    r <- M$engine$solve(M$model, data.frame(BW = bw), ev, tt, rtol = 1e-5, atol = 1e-10)
    w <- tt >= t_mdz
    list(r = r, t = tt, w = w, t_mdz = t_mdz,
         auc = auc_trap(tt[w] - t_mdz, r$MDZ[w, 1]), cmax = max(r$MDZ[w, 1]),
         cyp_l = r$CYP3A_LIVER[which(w)[1], 1], cyp_g = r$CYP3A_GUT[which(w)[1], 1])
  }

  params <- reactive({
    shiny::req(input$mdz_dose > 0, input$bw > 0)
    list(rif = input$rif_dose, days = input$rif_days, route = input$mdz_route, dose = input$mdz_dose,
         washout = input$washout, bw = input$bw)
  }) |> debounce(500)

  sim <- reactive({
    p <- params()
    withProgress(message = "Simulating rifampicin and midazolam", value = 0.3, {
      ddi <- run(p$rif, p$days, p$route, p$dose, p$washout, p$bw)
      ctl <- run(0, 1, p$route, p$dose, 0, p$bw)
    })
    list(ddi = ddi, ctl = ctl, p = p)
  })

  output$engine_note <- renderText(sprintf("  (%s)", M$engine$name))

  output$aucr <- renderValueBox({
    s <- sim(); v <- s$ddi$auc / s$ctl$auc
    stat_box(sprintf("%.2f", v), "", "Midazolam AUC ratio", sprintf("%.0f%% change vs no rifampicin", 100 * (v - 1)), "chart-area",
             color = if (v < 0.2) "red" else if (v < 0.5) "yellow" else "blue")
  })
  output$cmaxr <- renderValueBox({
    s <- sim(); stat_box(sprintf("%.2f", s$ddi$cmax / s$ctl$cmax), "", "Cmax ratio", "With vs without rifampicin", "arrow-up")
  })
  output$liver <- renderValueBox({
    s <- sim(); stat_box(sprintf("%.2f", s$ddi$cyp_l), "x", "Hepatic CYP3A", "When midazolam is given", "vial")
  })
  output$gut <- renderValueBox({
    s <- sim(); stat_box(sprintf("%.2f", s$ddi$cyp_g), "x", "Intestinal CYP3A", "When midazolam is given", "vial")
  })

  output$mdz_plot <- renderPlot({
    s <- sim()
    lv <- c("With rifampicin", "Without")
    d <- rbind(data.frame(time = s$ddi$t[s$ddi$w] - s$ddi$t_mdz, value = s$ddi$r$MDZ[s$ddi$w, 1], series = lv[1]),
               data.frame(time = s$ctl$t[s$ctl$w] - s$ctl$t_mdz, value = s$ctl$r$MDZ[s$ctl$w, 1], series = lv[2]))
    if (isTRUE(input$log)) d$value <- ifelse(d$value > max(d$value) * 1e-4, d$value, NA)
    p <- series_plot(d, "Hours after the midazolam dose", "Midazolam (ng/mL)", lv)
    if (isTRUE(input$log)) p + scale_y_log10(labels = plain_number) else p
  })

  output$cyp_plot <- renderPlot({
    s <- sim(); r <- s$ddi$r
    lv <- c("Liver", "Gut wall")
    d <- rbind(data.frame(time = s$ddi$t / 24, value = r$CYP3A_LIVER[, 1], series = lv[1]),
               data.frame(time = s$ddi$t / 24, value = r$CYP3A_GUT[, 1], series = lv[2]))
    series_plot(d, "Day", "CYP3A activity (x baseline)", lv) +
      geom_vline(xintercept = s$ddi$t_mdz / 24, colour = PAL$ink_3, linetype = "22") +
      geom_hline(yintercept = 1, colour = PAL$ink_3) +
      labs(subtitle = "Dashed line: midazolam dose")
  })

  output$rif_plot <- renderPlot({
    s <- sim()
    ggplot(data.frame(time = s$ddi$t / 24, v = s$ddi$r$RIF[, 1]), aes(time, v)) +
      geom_line(colour = PAL$blue_ink, linewidth = 1) +
      labs(x = "Day", y = "Rifampicin (ug/mL)") + theme_sim(12)
  })

  output$tbl <- DT::renderDataTable({
    s <- sim()
    small_table(data.frame(
      Metric = c("Midazolam AUC 0-24 h (ng*h/mL)", "Midazolam Cmax (ng/mL)", "AUC ratio", "Cmax ratio"),
      `Without rifampicin` = c(sprintf("%.1f", s$ctl$auc), sprintf("%.2f", s$ctl$cmax), "1", "1"),
      `With rifampicin` = c(sprintf("%.1f", s$ddi$auc), sprintf("%.2f", s$ddi$cmax),
                            sprintf("%.3f", s$ddi$auc / s$ctl$auc), sprintf("%.3f", s$ddi$cmax / s$ctl$cmax)),
      check.names = FALSE))
  })

  dr <- eventReactive(input$run_dr, {
    p <- params()
    doses <- c(5, 10, 25, 75, 150, 300, 600)
    ctl <- run(0, 1, p$route, p$dose, 0, p$bw)
    withProgress(message = "Rifampicin dose-response", value = 0, {
      do.call(rbind, lapply(doses, function(dd) {
        incProgress(1 / length(doses), detail = sprintf("%d mg", dd))
        x <- run(dd, p$days, p$route, p$dose, p$washout, p$bw)
        data.frame(dose = dd, aucr = x$auc / ctl$auc, liver = x$cyp_l, gut = x$cyp_g)
      }))
    })
  })

  output$dr_plot <- renderPlot({
    d <- dr()
    ggplot(d, aes(dose, aucr)) +
      geom_line(colour = PAL$blue_ink, linewidth = 1.1) + geom_point(colour = PAL$blue_ink, size = 3) +
      scale_x_log10(breaks = d$dose) + scale_y_continuous(limits = c(0, 1)) +
      labs(x = "Rifampicin daily dose (mg, log scale)", y = "Midazolam AUC ratio") + theme_sim(12)
  })

  output$dr_tbl <- DT::renderDataTable({
    d <- dr()
    small_table(data.frame(`Rifampicin (mg/day)` = d$dose, `AUC ratio` = sprintf("%.3f", d$aucr),
                           `Hepatic CYP3A (x)` = sprintf("%.2f", d$liver), `Gut CYP3A (x)` = sprintf("%.2f", d$gut),
                           check.names = FALSE))
  })
}

shinyApp(ui, server)
