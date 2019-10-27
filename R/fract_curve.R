#' Finding a point of fracture on a Kaplan-Meier curve.
#'
#' @import survival
#' @import survminer
#' @import tidyr
#' @import dplyr
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#' @importFrom graphics plot
#' @importFrom stats as.formula
#' @importFrom stats resid
#' @importFrom stats time
#' @importFrom utils write.csv
#'
#' @param data A data.frame-class object.
#' @param var.time A column name specifies the time variable in the data.
#' @param var.event A column name specifies the event variable in the data.
#' @param km.fit As an alternative to "data", "var.time" and "var.event", you can choose input a survfit-class object.
#' @param dir.output  The directory to output results.
#' @param get.df_of_IYs A logical if IYs will be obtained as data.frame in result.
#' @param fn.plot_pdf   The pdf file name for outputs.
#' @param fn.df_of_IYs The csv file name for IY output.
#'
#' @export

fract_curve <- function(
  data,
  var.time,
  var.event,
  km.fit=NULL,
  dir.output=NULL,
  fn.output=NULL,
  get.df_of_IYs =TRUE,
  fn.df_of_IYs=NULL
  ){

  if(get.df_of_IYs){
    if(is.null(fn.df_of_IYs)){stop(
      "Set 'FALSE' get.df_of_IYs argument or specify a name of csv to output."
    )
      }
    }


  if(is.null(km.fit)){

    if(
      is.null(var.time)|is.null(var.event)
    ){
      stop('"var.time" or "var.event" is not specified.')
      }

    data_surv <- data
    data_surv$dummy <- 1

    data_surv$val.time <- data_surv[ ,var.time]
    data_surv$val.event <- data_surv[ ,var.event]

    print(data_surv$val.time); print(data_surv$val.event)

    data_surv$Surv.obj <-
      with(
        data_surv,
        Surv(
          time  = val.time,
          event = val.event
          )
      )

    tmp_formula <- as.formula(
      sprintf(
        "%s ~ %s",
        "Surv.obj",
        paste(var.x,  collapse = " + ")
        )
      )

    km.fit <- survminer::surv_fit(
      # survminer::surv_fit()
      #   https://github.com/kassambara/survminer/issues/283
      tmp_formula,
      data= data_surv
      )
    }

  df.km.fit <- data.frame(
    time=km.fit$time,
    surv=km.fit$surv
    )

  # Obtain y-intersects of the secant lines -----------
  # on the scaled-KM curve.

  df.I.Y <- data.frame(
    km.fit$time,
    km.fit$surv
    ) %>%
    dplyr::mutate(
      I.Y =
        (km.fit.time-min(km.fit.time))/
        (max(km.fit.time)-min(km.fit.time)) +
        (km.fit.surv-min(km.fit.surv))/
        (max(km.fit.surv)-min(km.fit.surv)),
      rank.I.Y = rank(-I.Y)
      )

  # make secant line of the curve -----------------------
  # intersects the points of (0,1) and at which the value a is maximized..

  # linear regression between two points
  # returns a line intersects these two points.

  fit.IY <- stats::lm(
    I.Y ~ km.fit.time,
    df.I.Y %>%
      dplyr::filter(
        km.fit.time ==
          df.I.Y[
            df.I.Y$I.Y==max(df.I.Y$I.Y),
            "km.fit.time"
            ] |
          km.fit.time == min(km.fit.time)
        )
    )


  # Difference between the scaled-KM curve and the secant line.---------

  df.I.Y$pred.lm  <-
    stats::predict(fit.IY, df.I.Y)

  df.I.Y$resid <-
    df.I.Y$pred.lm - df.I.Y$I.Y


  # make secant line of the curve -----------------------
  # intersects the points of (0,1) and at which the value a is maximized..

  # linear regression between two points
  # returns a line intersects these two points.


  fit.km.fit.time <- stats::lm(
    surv ~ time,
    df.km.fit %>%
      dplyr::filter(
        time ==
          df.I.Y[
            df.I.Y$I.Y==max(df.I.Y$I.Y),
            "km.fit.time"
            ] |
          time==min(time)
      )
  )

  df.km.fit$pred.lm  <-
    stats::predict(fit.km.fit.time, df.km.fit)

  df.km.fit$resid <-
    df.km.fit$pred.lm - df.km.fit$surv



  ggdata <- df.I.Y %>%
    ggplot2::ggplot(
      ggplot2::aes(
        x=km.fit.time,
        y=I.Y
        )
      )

  ggdata_resid <- df.I.Y %>%
    ggplot2::ggplot(
      ggplot2::aes(
        x = km.fit.time,
        y=resid
      )
    )

  ggdata.km.fit_resid <- df.km.fit %>%
    ggplot2::ggplot(
      ggplot2::aes(
        x=time,
        y=resid
        )
    )



# Output ------------------------------------------------------------------

  # CSV filr of the Y-intercepts of secant lines.

  if(get.df_of_IYs){
    write.csv(
      sprintf(
        '%s/%s',
        dir.output,
        fn.df_of_IYs
        )
      )
    }

  # PDF file.

  pdf(
    sprintf(
      "%s/%s.pdf",
      dir.output,
      fn.plot_pdf
      ),
    height = 7,
    width = 10
  )
  plot(
    ggdata +
      ggplot2::geom_point() +
      ggplot2::geom_vline(
        xintercept =  df.I.Y[
          df.I.Y$I.Y==max(df.I.Y$I.Y),
          "km.fit.time"
          ]
      ) +
      scale_x_continuous(
        breaks = c(
          unname(
            stats::quantile(df.I.Y$km.fit.time,c(0, 0.25, 0.75, 1))
            ),
          df.I.Y[
            df.I.Y$I.Y==max(df.I.Y$I.Y),
            "km.fit.time"
            ]
        )
      ) +
      ggplot2::theme_bw()
  )
  plot(
    ggdata.km.fit_resid +
      ggplot2::geom_point() +
      ggplot2::geom_vline(
        xintercept =  df.I.Y[
          df.I.Y$I.Y==max(df.I.Y$I.Y),
          "km.fit.time"
          ]
      ) +
      scale_x_continuous(
        breaks = c(
          unname(stats::quantile(df.I.Y$km.fit.time,c(0, 0.25, 0.75, 1))),
          df.I.Y[
            df.I.Y$I.Y==max(df.I.Y$I.Y),
            "km.fit.time"
            ]
        )
      ) +
      ggplot2::theme_bw()
  )
  plot(
    ggdata_resid +
      ggplot2::geom_point() +
      ggplot2::theme_bw()
    )
  dev.off()
  }
# End runt---------



