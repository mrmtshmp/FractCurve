#' Decide the number of clusters by the fracture point of the height.
#'
#' @import survival
#' @import survminer
#' @import tidyr
#' @import dplyr
#' @import ggplot2
#' @importFrom stats hclust
#' @importFrom stats dist
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#' @importFrom graphics plot
#' @importFrom stats as.formula
#' @importFrom stats resid
#' @importFrom stats time
#' @importFrom utils write.csv
#' @importFrom tidyr expand_grid
#' @importFrom vegan vegdist
#' @importFrom coda.base dist
#' @importFrom NbClust NbClust
#'
#' @param df.features A data.frame-class object.
#' @param df.phenotype A data.frame-class object.
#' @param method.dist.row A distance metric for clustering of features.
#' @param method.dist.col A distance metric for clustering of samples.
#' @param method.hclust.row An aggromerizing algorithm for clustering of features.
#' @param method.hclust.col An aggromerizing algorithm for clustering of samples.
#' @param dir.output  The directory to output results.
#' @param get.df_of_IYs A logical if IYs will be obtained as data.frame in result.
#' @param fn.plot_pdf   The pdf file name for outputs.
#' @param fn.df_of_IYs The csv file name for IY output.
#'
#' @export

fract_curve_clust <- function(
  df.features  = OTUdata_all,
  df.phenotype = pData(obj.ADS),
  var.phenoGroup = 'Disease',
  method.dist.row   ='manhattan',
  method.dist.col   ='manhattan',
  method.hclust.row ='ward',
  method.hclust.col ='ward',
  fisher_test=TRUE,
  dir.output = NULL,
  get.df_of_IYs,
  fn.plot_pdf = NULL,
  fn.df_of_IYs = NULL
  ){


  if(method.dist.col == "aitchison"){my.dist.col <- coda.base::dist}else{
    my.dist.col <- vegan::vegdist
  }

  if(method.dist.row == "aitchison"){my.dist.row <- coda.base::dist}else{
    my.dist.row <- vegan::vegdist
  }


  rowv.hc <- hclust(
    my.dist.row(
      as.matrix(
        df.features
        ),
      method=method.dist.row
      ),
    method = method.hclust.row
    )

  rowv <- rowv.hc %>%
    as.dendrogram()


  colv.hc <- hclust(
    my.dist.col(
      t(
        as.matrix(
          df.features
          )
        ),
      method = method.dist.col
      ),
    method = method.hclust.col
    )

  colv <- colv.hc %>%
    as.dendrogram()

  res.clust <- list(colv.hc, rowv.hc)

  df.height <-
    data.frame(
      height=max(res.clust[[1]]$height)-res.clust[[1]]$height,
      event =1
      )

  df.res.fracrcurve <- FractCurve::fract_curve(
    df.height,
    var.time = 'height',
    var.event = 'event',
    fn.plot_pdf = fn.plot_pdf,
    dir.output = dir.output,
    get.df_of_IYs = get.df_of_IYs,
    fn.df_of_IYs = fn.df_of_IYs
    )

  df.res.fracrcurve <-
    df.res.fracrcurve[
      order(df.res.fracrcurve$km.fit.time),
      ]

  clusters <- cutree(
    res.clust[[1]],
    k = which(df.res.fracrcurve$rank.I.Y==1)
    )

  if(fisher_test){

    group <- df.phenotype[,var.phenoGroup]

    names(group) <-
      rownames(
        df.phenotype
      )

    row.select <-
      expand_grid(
        a=unique(clusters),
        b=unique(clusters)
      ) %>%
      filter(a < b)

    res.test <- list()

    for(i in 1:nrow(row.select)){

      row.select_i <- unlist(c(row.select[i,'a'], row.select[i,'b']))

      print(clusters)
      print(clusters)

      res.test_i <- fisher.test(
        as.matrix(
          table(clusters, group)[
            row.select_i,
            ]
        )
      )
      res.test_i$vs.info <- sprintf("%s vs %s", row.select[i,'a'], row.select[i,'b'])
      res.test[[i]] <- res.test_i
    }

    df.res.test <- sapply(
      lapply(
        res.test,
        function(x){t(data.frame(unlist(x)))}
      ),
      function(x){unlist(data.frame(x))}
    ) %>%
      t() %>%
      data.frame() %>%
      mutate(
        p.value=round(as.numeric(p.value), 3),
        conf.int1=round(as.numeric(conf.int1), 3),
        conf.int2=round(as.numeric(conf.int2), 3),
        estimate.odds.ratio=round(as.numeric(estimate.odds.ratio), 3)
      )
    return(
      list(
        df.res.fracrcurve,
        res.clust,
        clusters,
        df.res.test
      )
    )
  }else{
    return(
      list(
        df.res.fracrcurve,
        res.clust,
        clusters
      )
    )
    }
  }


