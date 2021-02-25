# #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
#
#                      Copyright (c) 2018, 2019, 2020
#            Marcel Ribeiro-Dantas <marcel.ribeiro-dantas@curie.fr>
#
# This script is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# This script is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #


# GSE17618 ----------------------------------------------------------------

# Loading variables from previous CPU-intensive tasks ---------------------

setwd('GSE17618/')
load('GSE17618.RData')


# Functions for survival analysis and plotting ----------------------------

# GSE17618 has only 44 samples from patients (cell lines have obviously no
# survival data). Instead of inferring a new network, we used the regulatory
# network inferred from the largest cohort we have, that is, GSE34620 with 117
# patients. This is referred here by the variable RegScores_inserted_in_skeleton
# which is using the skeleton from the regulatory network inferred from GSE34620

# This function outputs a CSV file with information regarding the provided
# regulon for every one of the 44 patients. That includes survival and
# regulatory information.
survival_regulon <- function(tf_name, log=FALSE) {
  if (log==FALSE) {
    regulon <- as.data.frame(RegScores_inserted_in_skeleton$dif[, tf_name])
  } else {
    regulon <- as.data.frame(RegScores_inserted_in_skeleton_log$dif[, tf_name])
  }
  regulon$genename <- tf_name
  colnames(regulon) <- c('expression', 'genename')

  regulon_final <- regulon
  regulon_final$sample <- rownames(regulon)
  regulon_final <- merge(x=amostras, y=regulon_final, by.x='sample')
  regulon_final <- cbind(regulon_final$sample, regulon_final$genename, regulon_final$expression, regulon_final$`OS (months)`, regulon_final$status)
  colnames(regulon_final) <- c('amostra', 'genename', 'expression', 'meses', 'status')
  regulon_final <- as.data.frame(regulon_final)
  if (log==FALSE) {
    write.csv(regulon_final, file = paste(tf_name, '.csv', sep=''))
  } else {
    write.csv(regulon_final, file = paste(tf_name, '_log.csv', sep=''))
  }

}

library(readr)
library(survival)
library(dplyr)

setwd('/tmp/')

ggsurv <- function(s, CI = 'def', plot.cens = T, surv.col = 'gg.def',
                   cens.col = 'red', lty.est = 1, lty.ci = 2,
                   cens.shape = 3, cens.size = 1, back.white = F, xlab = 'Time',
                   ylab = 'Survival', main = ''){

  library(ggplot2)
  strata <- ifelse(is.null(s$strata) ==T, 1, length(s$strata))
  stopifnot(length(surv.col) == 1 | length(surv.col) == strata)
  stopifnot(length(lty.est) == 1 | length(lty.est) == strata)

  ggsurv.s <- function(s, CI = 'def', plot.cens = T, surv.col = 'gg.def',
                       cens.col = 'red', lty.est = 1, lty.ci = 2,
                       cens.shape = 3, back.white = F, xlab = 'Time',
                       ylab = 'Survival', main = ''){

    dat <- data.frame(time = c(0, s$time),
                      surv = c(1, s$surv),
                      up = c(1, s$upper),
                      low = c(1, s$lower),
                      cens = c(0, s$n.censor))
    dat.cens <- subset(dat, cens != 0)

    col <- ifelse(surv.col == 'gg.def', 'black', surv.col)

    pl <- ggplot(dat, aes(x = time, y = surv)) +
      xlab(xlab) + ylab(ylab) + ggtitle(main) +
      geom_step(col = col, lty = lty.est)

    pl <- if(CI == T | CI == 'def') {
      pl + geom_step(aes(y = up), color = col, lty = lty.ci) +
        geom_step(aes(y = low), color = col, lty = lty.ci)
    } else (pl)

    pl <- if(plot.cens == T & length(dat.cens) > 0){
      pl + geom_point(data = dat.cens, aes(y = surv), shape = cens.shape,
                      col = cens.col,
                      size = cens.size)
    } else if (plot.cens == T & length(dat.cens) == 0){
      stop ('There are no censored observations')
    } else(pl)

    pl <- if(back.white == T) {pl + theme_bw()
    } else (pl)
    pl
  }

  ggsurv.m <- function(s, CI = 'def', plot.cens = T, surv.col = 'gg.def',
                       cens.col = 'red', lty.est = 1, lty.ci = 2,
                       cens.shape = 3, back.white = F, xlab = 'Time',
                       ylab = 'Survival', main = '') {
    n <- s$strata

    groups <- factor(unlist(strsplit(names
                                     (s$strata), '='))[seq(2, 2*strata, by = 2)])
    gr.name <-  unlist(strsplit(names(s$strata), '='))[1]
    gr.df <- vector('list', strata)
    ind <- vector('list', strata)
    n.ind <- c(0,n); n.ind <- cumsum(n.ind)
    for(i in 1:strata) ind[[i]] <- (n.ind[i]+1):n.ind[i+1]

    for(i in 1:strata){
      gr.df[[i]] <- data.frame(
        time = c(0, s$time[ ind[[i]] ]),
        surv = c(1, s$surv[ ind[[i]] ]),
        up = c(1, s$upper[ ind[[i]] ]),
        low = c(1, s$lower[ind[[i]]]),
        cens = c(0, s$n.censor[ind[[i]]]),
        group = rep(groups[i], n[i] + 1))
    }

    dat <- do.call(rbind, gr.df)
    dat.cens <- subset(dat, cens != 0)

    pl <- ggplot(dat, aes(x = time, y = surv, group = group)) +
      xlab(xlab) + ylab(ylab) + ggtitle(main) +
      geom_step(aes(col = group, lty = group))

    col <- if (length(surv.col == 1)) {
      scale_colour_manual(name = gr.name, values = rep(surv.col, strata))
    } else{
      scale_colour_manual(name = gr.name, values = surv.col)
    }

    pl <- if (surv.col[1] != 'gg.def') {
      pl + col
    } else {
      pl + scale_colour_discrete(name = gr.name)
    }

    line <- if (length(lty.est) == 1) {
      scale_linetype_manual(name = gr.name, values = rep(lty.est, strata))
    } else {
      scale_linetype_manual(name = gr.name, values = lty.est)
    }

    pl <- pl + line

    pl <- if (CI == T) {
      if (length(surv.col) > 1 && length(lty.est) > 1) {
        stop(
          'Either surv.col or lty.est should be of length 1 in order
          to plot 95% CI with multiple strata'
        )
      } else if ((length(surv.col) > 1 | surv.col == 'gg.def')[1]) {
        pl + geom_step(aes(y = up, color = group), lty = lty.ci) +
          geom_step(aes(y = low, color = group), lty = lty.ci)
      } else{
        pl +  geom_step(aes(y = up, lty = group), col = surv.col) +
          geom_step(aes(y = low, lty = group), col = surv.col)
      }
    } else {
      pl
    }


    pl <- if (plot.cens == T & length(dat.cens) > 0) {
      pl + geom_point(data = dat.cens,
                      aes(y = surv, col = status),
                      shape = cens.shape,
                      col = cens.col,
                      size = cens.size)
    } else if (plot.cens == T & length(dat.cens) == 0) {
      stop ('There are no censored observations')
    } else
      (pl)

    pl <- if (back.white == T) {
      pl + theme_bw()
    } else
      (pl)
    pl
  }
  pl <- if (strata == 1) {
    ggsurv.s(
      s,
      CI ,
      plot.cens,
      surv.col ,
      cens.col,
      lty.est,
      lty.ci,
      cens.shape,
      back.white,
      xlab,
      ylab,
      main
    )
  } else {
    ggsurv.m(
      s,
      CI,
      plot.cens,
      surv.col ,
      cens.col,
      lty.est,
      lty.ci,
      cens.shape,
      back.white,
      xlab,
      ylab,
      main
    )
  }
  pl
}

.internal_format <- function(data) {
  pvalue <- pchisq(data$chisq, length(data$n)-1, lower.tail = FALSE)
  row <- data.frame('Down' = data$n[1], 'DownDecease' = data$obs[1], 'Up' = data$n[2],
                    'UpDecease' = data$obs[2], 'pvalue_logrank' = pvalue)
}

surv_log_rank <- function(.data) {
  log_result <- .data %>% do(
    survival::survdiff(formula = survival::Surv(os_months, os_status) ~ groupname, data = .) %>% .internal_format()
  )
  log_result
}

plot_km <- function(genename) {
  gene_target <- genename
  gene_data <- sarcoma %>% filter(genename == gene_target)
  gene_fit <- survfit(Surv(os_months, os_status) ~ groupname, data = gene_data)

  gene_plot <- ggsurv(gene_fit, xlab = 'Months Survival', ylab = 'Surviving',
                      main = paste(gene_target), cens.col = 'black', cens.size = 0.5) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.05), limits = c(0, 1)) + theme_classic() +
    theme(legend.title = element_blank(), legend.position = c(0.8, 0.95), legend.direction = 'horizontal' )

  gene_plot
}

# Running analysis ---------------------------------------------------------------------

# Setting output directory
setwd('/home/mribeirodantas/Dropbox/Paper Mestrado/Submission 2020/Final/scripts/Análise de Sobrevida/GSE17618/')

# Output files with regulon/survival information for all patients separated by
# regulator.
survival_regulon('ARNT2')
survival_regulon('PAX7')
survival_regulon('RUNX3')
survival_regulon('CREB3L1')
survival_regulon('MEF2C')
survival_regulon('PBX3')
survival_regulon('GLI3')

# Merge all CSV files for all the 7 putative master regulators
sarcoma <- list.files(pattern = "\\.csv$", full.names = TRUE) %>%
  lapply(read_csv) %>%
  bind_rows

sarcoma$X1 <- NULL

colnames(sarcoma) <- c('sample', 'genename', 'expression', 'os_months', 'os_status')
unique(sarcoma$genename)

# It's Up if the expression is bigger or equal the group median. down otherwise.
sarcoma <- sarcoma %>% group_by(genename) %>%
  mutate(
    group_median = median(expression),
    groupname = ifelse(expression <= group_median, 'Down', 'Up'),
    os_status = ifelse(os_status == "Dead", 1, 0)
  )

# Remove variable
sarcoma$group_median <- NULL

testable <- table(sarcoma[, c("genename", "groupname")]) %>% as.data.frame()
testable <- data.table::dcast(testable, genename ~ groupname, value.var = "Freq")
testable <- testable %>% filter(Up > 0, Down > 0)

survival <- sarcoma %>% filter(genename %in% testable$genename) %>%
  group_by(genename) %>% surv_log_rank()

survival$padjust <- p.adjust(p=survival$pvalue_logrank, method="fdr", n=length(survival$pvalue_logrank))
survival[survival$pvalue_logrank<0.05, ]

# Output survival information
# write.csv(survival, file = 'survival.csv')

# Plotting KMs ------------------------------------------------------------
png(filename = 'PNGs/ARNT2', width = 2000, height = 800, res = 200)
plot_km('ARNT2')
dev.off()
png(filename = 'PNGs/PBX3', width = 2000, height = 800, res = 200)
plot_km('PBX3')
dev.off()
png(filename = 'PNGs/MEF2C', width = 2000, height = 800, res = 200)
plot_km('MEF2C')
dev.off()
png(filename = 'PNGs/GLI3', width = 2000, height = 800, res = 200)
plot_km('GLI3')
dev.off()
png(filename = 'PNGs/RUNX3', width = 2000, height = 800, res = 200)
plot_km('RUNX3')
dev.off()
png(filename = 'PNGs/PAX7', width = 2000, height = 800, res = 200)
plot_km('PAX7')
dev.off()
png(filename = 'PNGs/CREB3L1', width = 2000, height = 800, res = 200)
plot_km('CREB3L1')
dev.off()


# GSE63157 ----------------------------------------------------------------
setwd('../GSE63157/')

# Merge all CSV files for all the 7 putative master regulators
sarcoma <- list.files(pattern = "\\.csv$", full.names = TRUE) %>%
  lapply(read_csv) %>%
  bind_rows

sarcoma$X1 <- NULL

colnames(sarcoma) <- c('sample', 'genename', 'expression', 'os_months', 'os_status')
unique(sarcoma$genename)

# It's Up if the expression is bigger or equal the group median. down otherwise.
sarcoma <- sarcoma %>% group_by(genename) %>%
  mutate(
    group_median = median(expression),
    groupname = ifelse(expression <= group_median, 'Down', 'Up'),
    os_status = ifelse(os_status == "Dead", 1, 0)
  )

# Remove variable
sarcoma$group_median <- NULL

testable <- table(sarcoma[, c("genename", "groupname")]) %>% as.data.frame()
testable <- data.table::dcast(testable, genename ~ groupname, value.var = "Freq")
testable <- testable %>% filter(Up > 0, Down > 0)

survival <- sarcoma %>% filter(genename %in% testable$genename) %>%
  group_by(genename) %>% surv_log_rank()

survival$padjust <- p.adjust(p=survival$pvalue_logrank, method="fdr", n=length(survival$pvalue_logrank))
survival[survival$pvalue_logrank<0.05, ]

# Output survival information
# write.csv(survival, file = 'survival.csv')

# Plotting KMs ------------------------------------------------------------
png(filename = 'PNGs/ARNT2', width = 2000, height = 800, res = 200)
plot_km('ARNT2')
dev.off()
png(filename = 'PNGs/PBX3', width = 2000, height = 800, res = 200)
plot_km('PBX3')
dev.off()
png(filename = 'PNGs/MEF2C', width = 2000, height = 800, res = 200)
plot_km('MEF2C')
dev.off()
png(filename = 'PNGs/GLI3', width = 2000, height = 800, res = 200)
plot_km('GLI3')
dev.off()
png(filename = 'PNGs/RUNX3', width = 2000, height = 800, res = 200)
plot_km('RUNX3')
dev.off()
png(filename = 'PNGs/PAX7', width = 2000, height = 800, res = 200)
plot_km('PAX7')
dev.off()
png(filename = 'PNGs/CREB3L1', width = 2000, height = 800, res = 200)
plot_km('CREB3L1')
dev.off()
