# ============================================================================
# CAR-T cells: whole-body PBPK with tumour binding, expansion and killing
# (Singh et al. 2020). Preclinical: mouse physiology, mouse xenografts.
#
# Where do infused CAR-T cells go, how many reach the tumour, how many
# CAR-target complexes form per tumour cell, and is that enough to clear it?
# The app varies the dose, the tumour, the CAR's affinity and the antigen
# density, which the paper ties together.
# ============================================================================

library(shiny)
library(shinydashboard)
library(DT)
library(ggplot2)
library(scales)

shared <- function(f) source(if (file.exists(file.path("R", f))) file.path("R", f) else file.path("..", "..", "shared", f))
for (f in c("theme.R", "ode_engine.R", "app_helpers.R")) shared(f)

M <- load_model("cart_pbpk_pd_singh2020.cpp")

V_BLOOD <- 0.944          # mL, 25 g mouse
PER_ML_PER_NM <- 6.022e11 # molecules per mL in a 1 nM solution
KON <- 1.08e-12           # 1/((number/mL) h), fitted

# Organs, with the volumes the model uses (extravascular space)
ORGANS <- data.frame(
  name = c("Lung", "Spleen", "Liver", "Kidney", "Brain", "Gut", "Tumour", "Rest of body"),
  cmt = c("C_E_Lung", "C_E_Spleen", "C_E_Liver", "C_E_Kidney", "C_E_Brain", "C_E_GI", NA, "C_E_Other"),
  stringsAsFactors = FALSE
)

# The CAR-T products of Table 1, with the antigen density of the xenograft used
PRODUCTS <- list(
  bcma = list(label = "Anti-BCMA (RPMI-8226 myeloma)", kd = 10, taa = 12590, car = 5000, vol = 0.05),
  cd19 = list(label = "Anti-CD19 (CD19+ HeLa)", kd = 5, taa = 50000, car = 5000, vol = 0.05),
  egfr = list(label = "Anti-EGFR (U87 glioma)", kd = 40, taa = 30899, car = 5000, vol = 0.05)
)

ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = tags$div(tags$span("PBPK", style = "font-weight: bold; font-size: 24px;"),
                                   tags$span(" CAR-T biodistribution", style = "font-size: 18px;")),
                  titleWidth = 350),
  dashboardSidebar(
    width = 320, custom_css,
    sidebarMenu(
      menuItem("Tumour and dose", tabName = "main", icon = icon("bullseye")),
      menuItem("Where the cells go", tabName = "biod", icon = icon("sitemap")),
      menuItem("Affinity and antigen", tabName = "grid", icon = icon("th")),
      menuItem("About", tabName = "about", icon = icon("info-circle"))
    ),
    tags$div(
      style = "padding: 18px; padding-top: 8px;",
      tags$div(style = "background: #fdf3e3; border-radius: 6px; padding: 8px 10px; margin-bottom: 12px;",
               tags$strong("Preclinical model. "), "Mouse physiology and mouse xenograft data; not scaled to humans.",
               style = "font-size: 12px;"),
      tags$h4("CAR-T product"),
      selectInput("product", NULL, width = "100%",
                  choices = stats::setNames(names(PRODUCTS), vapply(PRODUCTS, `[[`, "", "label"))),
      fluidRow(
        column(6, numericInput("kd", "Affinity Kd (nM)", 10, min = 0.05, max = 5000, step = 1)),
        column(6, numericInput("car", "CAR per T cell", 5000, min = 100, max = 50000, step = 500))
      ),
      tags$h4("Tumour"),
      fluidRow(
        column(6, numericInput("taa", "Antigen per cell", 12590, min = 100, max = 1e6, step = 1000)),
        column(6, numericInput("vol", "Volume at dosing (mm3)", 50, min = 5, max = 500, step = 10))
      ),
      tags$h4("Dose"),
      sliderInput("dose", "Million CAR-T cells, IV", min = 0.1, max = 20, value = 10, step = 0.1),
      checkboxInput("compare", "Show the untreated tumour", TRUE),
      numericInput("days", "Days simulated", 28, min = 7, max = 56, step = 7, width = "50%")
    )
  ),
  dashboardBody(
    custom_css,
    tabItems(
      tabItem(
        tabName = "main",
        fluidRow(valueBoxOutput("tgi", width = 3), valueBoxOutput("cplx", width = 3),
                 valueBoxOutput("intumour", width = 3), valueBoxOutput("expand", width = 3)),
        fluidRow(
          box(title = tags$div(tags$strong("Tumour volume"),
                               tags$span(textOutput("engine_note", inline = TRUE),
                                         style = "color: #6b7078; font-size: 13px; font-weight: normal;")),
              status = "primary", solidHeader = TRUE, width = 7, plotOutput("tumour_plot", height = "360px"),
              muted("Killing acts through four transit compartments, so the tumour keeps growing for a few days after dosing before it turns over.")),
          box(title = "CAR-T cells", status = "primary", solidHeader = TRUE, width = 5,
              plotOutput("cells_plot", height = "360px"),
              muted("Blood falls within hours as cells distribute; cells that reach the tumour and meet antigen expand there."))
        ),
        fluidRow(
          box(title = "Dose-response", status = "primary", solidHeader = TRUE, width = 12,
              plotOutput("dose_plot", height = "280px"),
              muted("Tumour volume at the end of the simulation against dose, with the untreated tumour for reference. Below a threshold dose the tumour escapes; above it, it is cleared."))
        )
      ),
      tabItem(
        tabName = "biod",
        fluidRow(
          box(title = "CAR-T cells in each organ (extravascular space)", status = "primary", solidHeader = TRUE,
              width = 12, plotOutput("biod_plot", height = "420px"),
              muted("Lung takes up most cells first, then spleen and liver; the liver also clears them. Transmigration rates were fitted to radiolabelled CAR-T biodistribution in mice."))
        ),
        fluidRow(
          box(title = "Peak cell concentration by organ", status = "primary", solidHeader = TRUE, width = 12,
              DTOutput("biod_table"))
        )
      ),
      tabItem(
        tabName = "grid",
        fluidRow(
          box(title = "Tumour outcome against CAR affinity and antigen density", status = "primary",
              solidHeader = TRUE, width = 12, plotOutput("grid_plot", height = "440px"),
              muted("Each cell is one simulation at the selected dose: tumour volume at the end as a percentage of the untreated tumour. Tighter binding (low Kd) and more antigen both give more complexes per tumour cell, but the response saturates - beyond a point, higher affinity adds nothing, the paper's main conclusion for CAR design."))
        )
      ),
      about_tab(
        "CAR-T cells: whole-body PBPK with tumour binding and killing",
        "A multiscale PK-PD model that follows CAR-T cells through the body of a tumour-bearing
         mouse, forms CAR-target complexes in the tumour, and lets those complexes drive both
         CAR-T expansion and tumour-cell killing. It was built from in-vitro killing, cytokine
         and proliferation assays, radiolabelled biodistribution studies and tumour-growth
         inhibition studies in xenograft models.",
        list("Blood flow through lung, spleen, gut, liver, kidney, brain, tumour and the rest of the body, each with a vascular and an extravascular space",
             "Cells cross into tissue at a first-order transmigration rate fitted per organ, return by lymph through a lymph-node compartment, and are lost in the liver",
             "In the tumour: CAR-target complex formation (kon, koff) from CAR density per T cell and antigen density per tumour cell",
             "Complexes per tumour cell drive CAR-T expansion (Emax/EC50) and tumour killing (Emax/IC50, through four transit compartments)",
             "Tumour grows exponentially; 1 mL of tumour is 1e8 cells",
             "Mouse (25 g) physiology throughout: organ volumes and flows, and parameters fitted to mouse data"),
        tags$span("Singh AP, Zheng X, Lin-Schmidt X, Chen W, Carpenter TJ, Zong A, Wang W, Heald DL.
                   Development of a quantitative relationship between CAR-affinity, antigen abundance,
                   tumor cell depletion and CAR-T cell expansion using a multiscale systems PK-PD model.
                   mAbs 2020;12:1688616. Generated from the MLXTRAN model code in its Supplementary
                   Material by tools/build_cart_pbpk.py; parameter values from that code and Table 2."),
        M$engine$name
      )
    )
  )
)

server <- function(input, output, session) {

  observeEvent(input$product, {
    p <- PRODUCTS[[input$product]]
    updateNumericInput(session, "kd", value = p$kd)
    updateNumericInput(session, "car", value = p$car)
    updateNumericInput(session, "taa", value = p$taa)
    updateNumericInput(session, "vol", value = p$vol * 1000)
  }, ignoreInit = TRUE)

  pars <- reactive({
    shiny::req(input$kd > 0, input$car > 0, input$taa > 0, input$vol > 0, input$days >= 7)
    data.frame(KOFF = KON * input$kd * PER_ML_PER_NM, DENSITY_CAR = input$car,
               DENSITY_TAA = input$taa, VTUMOR0 = input$vol / 1000)
  }) |> debounce(500)

  dose_ev <- function(million, id = NA) {
    data.frame(time = 0, cmt = "C_Blood", amt = million * 1e6 / V_BLOOD, ID = id)
  }

  sim <- reactive({
    P <- pars()
    tt <- seq(0, 24 * input$days, by = 2)
    withProgress(message = "Simulating CAR-T biodistribution", value = 0.3, {
      treated <- M$engine$solve(M$model, P, dose_ev(input$dose), tt, rtol = 1e-6, atol = 1e-6)
      control <- M$engine$solve(M$model, P, NULL, tt, rtol = 1e-6, atol = 1e-6)
    })
    list(t = tt / 24, r = treated, ctrl = control)
  })

  output$engine_note <- renderText(sprintf("  (%s)", M$engine$name))

  output$tgi <- renderValueBox({
    x <- sim(); n <- length(x$t)
    tgi <- 100 * (1 - x$r$TumorVolume[n, 1] / x$ctrl$TumorVolume[n, 1])
    stat_box(sprintf("%.0f", tgi), " %", "Tumour growth inhibition",
             sprintf("Day %g, vs untreated", max(x$t)), "bullseye")
  })
  output$cplx <- renderValueBox({
    x <- sim()
    stat_box(signif(max(x$r$CplxPT[, 1]), 3), "", "Complexes per tumour cell", "Peak", "link")
  })
  output$intumour <- renderValueBox({
    x <- sim()
    stat_box(label_number(scale_cut = cut_short_scale(), accuracy = 0.1)(max(x$r$CARTtumor[, 1])), " /mL",
             "CAR-T in the tumour", "Peak concentration", "crosshairs")
  })
  output$expand <- renderValueBox({
    x <- sim()
    f <- max(x$r$CARTtumor[, 1]) / max(x$r$CARTtumor[x$t <= 1, 1])
    stat_box(sprintf("%.1f", f), " x", "Expansion in the tumour", "Peak vs day 1", "chart-line")
  })

  output$tumour_plot <- renderPlot({
    x <- sim()
    lv <- c(sprintf("%g million CAR-T cells", input$dose), "Untreated")
    d <- data.frame(time = x$t, value = x$r$TumorVolume[, 1], series = lv[1])
    if (isTRUE(input$compare)) d <- rbind(d, data.frame(time = x$t, value = x$ctrl$TumorVolume[, 1], series = lv[2]))
    series_plot(d, "Days after infusion", "Tumour volume (mm3)", lv[seq_len(length(unique(d$series)))])
  })

  output$cells_plot <- renderPlot({
    x <- sim()
    lv <- c("Blood (per uL)", "Tumour (per mL)")
    d <- rbind(data.frame(time = x$t, value = pmax(x$r$CARTblood[, 1], 1e-2), series = lv[1]),
               data.frame(time = x$t, value = pmax(x$r$CARTtumor[, 1], 1e-2), series = lv[2]))
    series_plot(d, "Days after infusion", "CAR-T cells", lv) +
      scale_y_log10(labels = label_number(scale_cut = cut_short_scale(), drop0trailing = TRUE))
  })

  output$dose_plot <- renderPlot({
    P <- pars()
    doses <- c(0.1, 0.3, 1, 3, 5, 10, 20)
    tt <- seq(0, 24 * input$days, by = 6)
    ev <- do.call(rbind, lapply(seq_along(doses), function(i) dose_ev(doses[i], id = i)))
    withProgress(message = "Simulating the dose range", value = 0.4, {
      r <- M$engine$solve(M$model, P[rep(1, length(doses)), ], ev, tt, rtol = 1e-6, atol = 1e-6)
      ctrl <- M$engine$solve(M$model, P, NULL, tt, rtol = 1e-6, atol = 1e-6)
    })
    n <- length(tt)
    d <- data.frame(dose = doses, value = r$TumorVolume[n, ])
    ggplot(d, aes(dose, value)) +
      geom_hline(yintercept = ctrl$TumorVolume[n, 1], colour = PAL$ink_3, linetype = "22") +
      geom_line(colour = PAL$blue_ink, linewidth = 1) + geom_point(colour = PAL$blue_ink, size = 2.6) +
      scale_x_log10(breaks = doses) +
      labs(x = "Million CAR-T cells (log scale)", y = sprintf("Tumour volume, day %g (mm3)", max(tt) / 24),
           subtitle = "Dashed: untreated") + theme_sim(12)
  })

  output$biod_plot <- renderPlot({
    x <- sim()
    org <- ORGANS[!is.na(ORGANS$cmt), ]
    d <- do.call(rbind, lapply(seq_len(nrow(org)), function(i) {
      data.frame(time = x$t, value = pmax(x$r$states[, org$cmt[i], 1], 1e-2), series = org$name[i])
    }))
    d <- rbind(d, data.frame(time = x$t, value = pmax(x$r$CARTtumor[, 1], 1e-2), series = "Tumour"))
    lv <- c(org$name, "Tumour")
    series_plot(d, "Days after infusion", "CAR-T cells per mL of tissue", lv) +
      scale_y_log10(labels = label_number(scale_cut = cut_short_scale(), drop0trailing = TRUE)) +
      guides(colour = guide_legend(ncol = 2))
  })

  output$biod_table <- renderDT({
    x <- sim()
    org <- ORGANS[!is.na(ORGANS$cmt), ]
    peak <- vapply(org$cmt, function(v) max(x$r$states[, v, 1]), 0)
    tmax <- vapply(org$cmt, function(v) x$t[which.max(x$r$states[, v, 1])], 0)
    d <- data.frame(Organ = c(org$name, "Tumour"),
                    `Peak cells per mL` = signif(c(peak, max(x$r$CARTtumor[, 1])), 3),
                    `Day of peak` = round(c(tmax, x$t[which.max(x$r$CARTtumor[, 1])]), 1),
                    check.names = FALSE)
    small_table(d[order(-d$`Peak cells per mL`), ])
  })

  output$grid_plot <- renderPlot({
    P <- pars()
    kds <- c(0.1, 1, 10, 100, 1000)
    taas <- c(1e3, 1e4, 5e4, 2e5, 1e6)
    g <- expand.grid(kd = kds, taa = taas)
    Pg <- P[rep(1, nrow(g)), ]
    Pg$KOFF <- KON * g$kd * PER_ML_PER_NM
    Pg$DENSITY_TAA <- g$taa
    tt <- seq(0, 24 * input$days, by = 6)
    ev <- do.call(rbind, lapply(seq_len(nrow(g)), function(i) dose_ev(input$dose, id = i)))
    withProgress(message = "Simulating the affinity-antigen grid", detail = "25 simulations", value = 0.3, {
      r <- M$engine$solve(M$model, Pg, ev, tt, rtol = 1e-5, atol = 1e-5)
      ctrl <- M$engine$solve(M$model, P, NULL, tt, rtol = 1e-5, atol = 1e-5)
    })
    n <- length(tt)
    g$pct <- pmin(100, 100 * r$TumorVolume[n, ] / ctrl$TumorVolume[n, 1])
    g$kdf <- factor(g$kd, levels = kds, labels = c("0.1", "1", "10", "100", "1000"))
    g$taaf <- factor(g$taa, levels = taas, labels = c("1e3", "1e4", "5e4", "2e5", "1e6"))
    ggplot(g, aes(kdf, taaf, fill = pct)) +
      geom_tile(colour = "white", linewidth = 2) +
      geom_text(aes(label = sprintf("%.0f%%", pct)), colour = PAL$ink, size = 4) +
      scale_fill_gradient(low = "#dceaf8", high = "#2a78d6", name = "% of untreated", limits = c(0, 100)) +
      labs(x = "CAR affinity Kd (nM), tighter to the left",
           y = "Antigen per tumour cell",
           subtitle = sprintf("Tumour volume on day %g as a percentage of untreated, at %g million CAR-T cells",
                              max(tt) / 24, input$dose)) +
      theme_sim(12)
  })
}

shinyApp(ui, server)
