# ---------------------------------------------------------------------------
# Shiny app launchers ----
# ---------------------------------------------------------------------------

#' Launch Interactive Threshold Explorer
#'
#' Opens a Shiny application that lets you drag a threshold slider and
#' immediately see how the confusion matrix and all metric values change.
#' The FOM curve is displayed with the currently selected threshold marked.
#'
#' @param toc A `toc_data` object (from [toc_data()] or [toc_from_scores()]).
#' @param ... Extra arguments passed to [shiny::shinyApp()].
#' @return Invisibly, the Shiny app object (so it can be embedded in R Markdown).
#' @export
#' @examples
#' \dontrun{
#' set.seed(1)
#' td <- toc_from_scores(runif(500), runif(500) > 0.4)
#' launch_threshold_explorer(td)
#' }
launch_threshold_explorer <- function(toc, ...) {
  rlang::check_installed("shiny", reason = "to run the threshold explorer app")

  P  <- attr(toc, "n_positive")
  Q  <- attr(toc, "n_negative")
  N  <- attr(toc, "n_total")
  pr <- P / N

  ui <- shiny::fluidPage(
    shiny::tags$head(shiny::tags$style(
      "body { font-family: 'Helvetica Neue', Helvetica, sans-serif; }"
    )),
    shiny::titlePanel("Threshold Explorer"),
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        width = 3,
        shiny::h4("Settings"),
        shiny::sliderInput(
          "k", "Threshold rank (k = # predicted positives)",
          min = 0, max = N, value = P, step = 1
        ),
        shiny::hr(),
        shiny::h4("Confusion Matrix"),
        shiny::tableOutput("cm_table"),
        shiny::hr(),
        shiny::h4("Agreement Metrics"),
        shiny::tableOutput("agreement_table"),
        shiny::hr(),
        shiny::h4("Skill Metrics"),
        shiny::tableOutput("skill_table")
      ),
      shiny::mainPanel(
        width = 9,
        shiny::fluidRow(
          shiny::column(6, shiny::plotOutput("fom_plot",   height = "380px")),
          shiny::column(6, shiny::plotOutput("roc_plot",   height = "380px"))
        ),
        shiny::fluidRow(
          shiny::column(6, shiny::plotOutput("oa_plot",    height = "280px")),
          shiny::column(6, shiny::plotOutput("gss_plot",   height = "280px"))
        )
      )
    )
  )

  server <- function(input, output, session) {
    cm_reactive <- shiny::reactive({
      k   <- as.integer(round(input$k))
      idx <- which(toc$k == k)
      if (length(idx) != 1) return(NULL)
      row <- toc[idx, , drop = FALSE]
      confusion_matrix(row$hits[[1]], row$fa[[1]], row$misses[[1]], row$cr[[1]])
    })

    output$cm_table <- shiny::renderTable({
      cm <- cm_reactive()
      shiny::req(cm)
      data.frame(
        ` `      = c("Predicted +", "Predicted -"),
        `Obs +`  = c(cm$hits,   cm$misses),
        `Obs -`  = c(cm$fa,     cm$cr),
        check.names = FALSE
      )
    }, rownames = FALSE, digits = 0)

    output$agreement_table <- shiny::renderTable({
      cm <- cm_reactive()
      shiny::req(cm)
      data.frame(
        Metric = c("OA", "BA", "MCC", "Kappa", "F1"),
        Value  = round(c(oa(cm), ba(cm), mcc(cm), kappa_score(cm), f1_score(cm)), 4)
      )
    }, rownames = FALSE)

    output$skill_table <- shiny::renderTable({
      cm <- cm_reactive()
      shiny::req(cm)
      data.frame(
        Metric = c("FOM/CSI", "GSS", "PSS", "HSS"),
        Value  = round(c(fom(cm), gss(cm), pss(cm), hss(cm)), 4)
      )
    }, rownames = FALSE)

    make_fom_plot <- function(k_sel) {
      p <- plot_fom(toc) +
        ggplot2::geom_vline(
          xintercept = k_sel / N,
          linetype   = "dashed",
          colour     = "orange",
          linewidth  = 1
        ) +
        ggplot2::geom_point(
          data = tibble::tibble(
            x = k_sel / N,
            y = toc$fom_curve[toc$k == k_sel]
          ),
          ggplot2::aes(x = x, y = y),
          colour = "orange", size = 4, inherit.aes = FALSE
        )
      p
    }

    output$fom_plot <- shiny::renderPlot({ make_fom_plot(input$k) })
    output$roc_plot <- shiny::renderPlot({
      plot_roc(toc) +
        ggplot2::geom_point(
          data = tibble::tibble(
            x = toc$fpr[toc$k == input$k],
            y = toc$tpr[toc$k == input$k]
          ),
          ggplot2::aes(x = x, y = y),
          colour = "orange", size = 4, inherit.aes = FALSE
        )
    })
    output$oa_plot  <- shiny::renderPlot({
      mc <- toc_metric_curve(toc, oa, "OA")
      plot_metric_curve(mc, "OA") +
        ggplot2::geom_vline(xintercept = input$k / N,
                            linetype = "dashed", colour = "orange", linewidth = 1)
    })
    output$gss_plot <- shiny::renderPlot({
      mc <- toc_metric_curve(toc, gss, "GSS")
      plot_metric_curve(mc, "GSS") +
        ggplot2::geom_vline(xintercept = input$k / N,
                            linetype = "dashed", colour = "orange", linewidth = 1)
    })
  }

  app <- shiny::shinyApp(ui = ui, server = server, ...)
  shiny::runApp(app)
  invisible(app)
}

# ---------------------------------------------------------------------------
# Bounds explorer ----
# ---------------------------------------------------------------------------

#' Launch Interactive Bounds Explorer
#'
#' Opens a Shiny application that plots the theoretical bounds on AUC and AFOM
#' as a function of the Maximum FOM (MFOM), with a slider to adjust prevalence
#' and observe how the feasible region changes.
#'
#' @param ... Extra arguments passed to [shiny::shinyApp()].
#' @return Invisibly, the Shiny app object.
#' @export
#' @examples
#' \dontrun{
#' launch_bounds_explorer()
#' }
launch_bounds_explorer <- function(...) {
  rlang::check_installed("shiny", reason = "to run the bounds explorer app")

  ui <- shiny::fluidPage(
    shiny::titlePanel("Bounds Explorer: AUC and AFOM vs Maximum FOM"),
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        width = 3,
        shiny::sliderInput(
          "prev", "Prevalence (P/N)",
          min = 0.01, max = 0.5, value = 0.1, step = 0.01
        ),
        shiny::hr(),
        shiny::p(
          "The shaded region shows the feasible space for a binary classifier
           given its Maximum FOM (x-axis) and the specified prevalence."
        ),
        shiny::p(
          "Points outside the bounds correspond to classifiers with an
           impossible combination of AUC and MFOM."
        )
      ),
      shiny::mainPanel(
        width = 9,
        shiny::fluidRow(
          shiny::column(6, shiny::plotOutput("auc_plot",  height = "400px")),
          shiny::column(6, shiny::plotOutput("afom_plot", height = "400px"))
        )
      )
    )
  )

  server <- function(input, output, session) {
    # AUC bounds as function of MFOM (from thesis formulas)
    auc_upper_fn <- function(mfom_val, pr) {
      # If MFOM >= pr: upper AUC uses triangle formula
      P_val <- 1          # normalised (P=1, E=1/pr)
      E_val <- 1 / pr
      a     <- mfom_val * P_val
      if (mfom_val >= pr) {
        num <- 2 * P_val * (E_val - P_val) - (P_val^2 / a - P_val) * (P_val - a)
        den <- 2 * P_val * (E_val - P_val)
      } else {
        num <- 2 * P_val * (E_val - P_val) - (E_val - P_val) * (P_val - a)
        den <- 2 * P_val * (E_val - P_val)
      }
      num / den
    }
    auc_lower_fn <- function(mfom_val) (1 + mfom_val) / 2

    output$auc_plot <- shiny::renderPlot({
      pr  <- input$prev
      seq_mfom <- seq(pr, 1, length.out = 200)

      bounds <- tibble::tibble(
        MFOM  = seq_mfom * 100,
        Upper = pmin(sapply(seq_mfom, auc_upper_fn, pr = pr), 1) * 100,
        Lower = pmax(sapply(seq_mfom, auc_lower_fn), pr) * 100
      )
      bounds$Upper[1] <- bounds$Lower[1] <- 50

      ggplot2::ggplot(bounds, ggplot2::aes(x = MFOM)) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = Lower, ymax = Upper),
                             fill = "steelblue", alpha = 0.15) +
        ggplot2::geom_line(ggplot2::aes(y = Upper), colour = "grey20",
                           linewidth = 1.4) +
        ggplot2::geom_line(ggplot2::aes(y = Lower), colour = "grey20",
                           linewidth = 1.4) +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                             colour = "grey50") +
        ggplot2::scale_x_continuous(limits = c(0, 100), expand = c(0, 0)) +
        ggplot2::scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
        ggplot2::labs(
          x        = "Maximum FOM / CSI (%)",
          y        = "AUC (%)",
          title    = "Feasible Region: AUC vs MFOM",
          subtitle = sprintf("Prevalence = %.2f", pr)
        ) +
        ggplot2::coord_fixed() +
        binm_theme()
    })

    output$afom_plot <- shiny::renderPlot({
      pr <- input$prev
      seq_mfom <- seq(pr, 1, length.out = 200)

      # AFOM upper: analytical formula from the thesis
      afom_upper_fn <- function(mfom_val, pr) {
        P_val  <- 1
        E_val  <- 1 / pr
        x_A    <- mfom_val * P_val
        # find x_B numerically
        xb_fn  <- function(x) P_val / x - mfom_val - pr * x / ((1 - pr) * x + P_val)
        x_B    <- tryCatch(
          stats::uniroot(xb_fn, lower = P_val, upper = E_val)$root,
          error = function(e) E_val
        )
        h_fn   <- function(x) pmin(x, (P_val - x_A) / (x_B - x_A) * (x - x_B) + P_val, P_val)
        u_fn   <- function(x) pr * x / ((1 - pr) * x + P_val)
        f_fn   <- function(x) h_fn(x) / (x - h_fn(x) + P_val)
        num    <- stats::integrate(f_fn, lower = 0, upper = E_val)$value -
                  P_val * pr / 2
        tot    <- P_val * (0.5 - log(pr) - pr / 2)
        afom_v <- num / tot
        aufom_v <- (stats::integrate(u_fn, lower = 0, upper = E_val)$value -
                    P_val * pr / 2) / tot
        (afom_v - aufom_v) / (1 - aufom_v)
      }

      upper_vals <- tryCatch(
        sapply(seq_mfom, afom_upper_fn, pr = pr),
        error = function(e) rep(NA, length(seq_mfom))
      )

      bounds <- tibble::tibble(
        MFOM  = seq_mfom * 100,
        Upper = pmax(upper_vals, 0) * 100
      )

      ggplot2::ggplot(bounds, ggplot2::aes(x = MFOM)) +
        ggplot2::geom_line(ggplot2::aes(y = Upper), colour = "grey20",
                           linewidth = 1.4, na.rm = TRUE) +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                             colour = "grey50") +
        ggplot2::scale_x_continuous(limits = c(0, 100), expand = c(0, 0)) +
        ggplot2::scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
        ggplot2::labs(
          x        = "Maximum FOM / CSI (%)",
          y        = "AFOM Skill (%)",
          title    = "Upper Bound: AFOM Skill vs MFOM",
          subtitle = sprintf("Prevalence = %.2f", pr)
        ) +
        ggplot2::coord_fixed() +
        binm_theme()
    })
  }

  app <- shiny::shinyApp(ui = ui, server = server, ...)
  shiny::runApp(app)
  invisible(app)
}
