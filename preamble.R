suppressMessages(library(tidyverse))
suppressMessages(library(ggbeeswarm))
suppressMessages(library(vroom))
suppressMessages(library(lubridate))
suppressMessages(library(scales))
suppressMessages(library(patchwork))
suppressMessages(library(binom))

# Shortcut functions
read.tsv = function(file){read.csv(file, sep = '\t')}
write.tsv = function(df, name){write.table(df, name, row.names = F, quote = F, sep = '\t')}
nl = theme(legend.position='none') # no legend
eb = element_blank() # for using the theme()

headn = function(data, n = 3    ){
    # Show the top of a table and display the total number of rows
    print(nrow(data))
    return(head(data, n))
    }

replace_na = function(data, x, default){
	x <- enquo(x)
	mutate(data, !!x := ifelse(is.na(!!x), default, !!x))
}

# Extract center of bins from cut()
bin_center = function(x){
    return(
        sapply(x, function(name){mean(as.numeric(unlist(strsplit(gsub("\\(|\\]", "", name), ","))))})
        )
    }

# Axes
pretty_log_breaks = function(x) unique(10^(seq(floor(min(log10(x))) - 1, ceiling(max(log10(x))) + 1, 1)))

xlog = function(lower = NULL, upper = NULL, expand = c(0,0), ...){ return(scale_x_continuous(limits = c(lower, upper), trans = log10_trans(), labels = trans_format("log10", math_format(10^.x)), expand = expand, ...))}
ylog = function(lower = NULL, upper = NULL, expand = c(0,0), ...){ return(scale_y_continuous(limits = c(lower, upper), trans = log10_trans(), labels = trans_format("log10", math_format(10^.x)), expand = expand, ...))}
logticks = function(s = 'bl', lw = 0.2, ...) annotation_logticks(sides = s, linewidth = lw, ...)

xcont = function(lower = NULL, upper = NULL, expand = c(0, 0), ...){ return(scale_x_continuous(limits = c(lower, upper), expand = expand, ...))}
ycont = function(lower = NULL, upper = NULL, expand = c(0, 0), ...){ return(scale_y_continuous(limits = c(lower, upper), expand = expand, ...))}

# Titled axis text
tilt_text = function(angle = -45, hjust = 0.1, vjust = 0.1, size = 8) theme(axis.text.x = element_text(angle = angle, hjust = hjust, vjust = vjust, size = size)) 

### Jupyter tricks ###

# Plot size
size = function(h,w){
	plot_h <<- h
	plot_w <<- w
	options(repr.plot.width=w, repr.plot.height=h, repr.plot.res=300)
}

size(2.5, 3.5)

# Display more columns
maxcols = function(n){ options(repr.matrix.max.cols=n) }
maxcols(150)

# Extra margin (for long legends that get cut out)
plot_margin = function(right = 0, left = 0, top = 0, bottom = 0) theme(plot.margin = unit(c(top, right, bottom, left), "lines"))

# Simple width and height function with sane defaults
wide = function(){size(plot_h, 8)}
narrow = function(){size(plot_h, 3.5)}
tall = function(){size(5, plot_w)}
short = function(){size(2.5, plot_w)}
big = function(){size(4.5, 8.5)}

library(tools)

psave <- function(name) {
  # Extract folder name from the path
  folder <- dirname(name)
  
  # Check if "PNG" folder exists, if not create it
  pdf_folder <- file.path(folder, "PDF")
  if (!file.exists(pdf_folder)) {
    dir.create(pdf_folder)
  }
    
  # Extract plot name
  plot_name <- file_path_sans_ext(basename(name))
    
  # Save plot as PNG
  png_path <- file.path(folder, paste0(plot_name, ".png"))
  ggsave(png_path, width=plot_w, height=plot_h, dpi = 300, device = "png")
  
  # Save plot as PNG in the "PNG" folder
  pdf_path <- file.path(pdf_folder, paste0(plot_name, ".pdf"))
  ggsave(pdf_path, width=plot_w, height=plot_h, dpi = 300, device = cairo_pdf)
}

psave_only=function(name){
	ggsave(name, width=plot_w, height=plot_h, dpi = 300)
}

######################
### Bayesian stuff ###
######################

library(HDInterval) # highest density credible intervals

# Closed-form Beta-binomial posterior parameter estimates. k and n may/may not be vectors. (1,1): uniform prior. (0.5,0.5): Jeffreys prior.
binom_def_prior_a = 0.5
binom_def_prior_b = 0.5
def_conf = 0.99

binom_post_a = function(k, n, prior_a = binom_def_prior_a, prior_b = binom_def_prior_b, ...) prior_a + sum(k)
binom_post_b = function(k, n, prior_a = binom_def_prior_a, prior_b = binom_def_prior_b, ...) prior_b + sum(n - k)
binom_post_mean = function(k, n, ...) binom_post_a(k, n, ...) / (binom_post_a(k, n, ...) + binom_post_b(k, n, ...))
binom_hdi_low = function(k, n, conf = def_conf, ...) hdi(qbeta, conf, shape1 = binom_post_a(k, n, ...), shape2 = binom_post_b(k, n, ...))[1]
binom_hdi_high = function(k, n, conf = def_conf, ...) hdi(qbeta, conf, shape1 = binom_post_a(k, n, ...), shape2 = binom_post_b(k, n, ...))[2]

# Closed-form Gamma-Poisson
poisson_def_prior_a = 0.0
poisson_def_prior_b = 0.0

poisson_post_a = function(k, prior_a = poisson_def_prior_a, prior_b = poisson_def_prior_b, ...) prior_a + sum(k)
poisson_post_b = function(k, prior_a = poisson_def_prior_a, prior_b = poisson_def_prior_b, ...) prior_b + length(k)
poisson_post_mean = function(k, prior_a = poisson_def_prior_a, prior_b = poisson_def_prior_b, ...) poisson_post_a(k, ...) / poisson_post_b(k, ...)
poisson_hdi_low = function(k, conf = 0.95, ...) hdi(qgamma, conf, shape = poisson_post_a(k, ...), rate = poisson_post_b(k, ...))[1]
poisson_hdi_high = function(k, conf = 0.95, ...) hdi(qgamma, conf, shape = poisson_post_a(k, ...), rate = poisson_post_b(k, ...))[2]

# Normal-inverse-gamma updating. Mean ~ N(mu, sigma), Var ~ IG(alpha, beta)
norm_post_NIG = function(x, prior_mu, prior_sigma, conf = def_conf, prior_alpha = 1){

    # Calculate the posterior parameters
    n <- length(x)
    x_bar <- mean(x)
    s_sq <- var(x)
    
    # Uninformative prior for beta: sample variance
    prior_beta = (n-1) * s_sq / 2
    prior_lambda <- 1 / (prior_sigma^2) # precision

    post_mu <- (prior_lambda * prior_mu + n * x_bar) / (prior_lambda + n)
    post_lambda <- prior_lambda + n
    post_alpha <- prior_alpha + n/2
    post_beta <- prior_beta + (n * s_sq + (prior_lambda * n * (x_bar - prior_mu)^2) / (prior_lambda + n)) / 2
    
    post_sigma = sqrt(1 / post_lambda)
    
    # Calculate the parameters of the marginal distribution of mu (integrating away sigma)
    scale <- sqrt(post_beta / (post_lambda * post_alpha))

    # Get the credible interval
    conf_tail = (1 - conf) / 2
    post_mu_low <- qt(conf_tail, df = 2 * post_alpha) * scale + post_mu
    post_mu_high <- qt(1 - conf_tail, df = 2 * post_alpha) * scale + post_mu

    return(c(post_mu, post_mu_low, post_mu_high))
}

norm_post_mu = function(x, prior_mu, prior_sigma, prior_alpha = 1) norm_post_NIG(x, prior_mu, prior_sigma, prior_alpha = prior_alpha)[1]
norm_post_low = function(x, prior_mu, prior_sigma, prior_alpha = 1) norm_post_NIG(x, prior_mu, prior_sigma, prior_alpha = prior_alpha)[2]
norm_post_high = function(x, prior_mu, prior_sigma, prior_alpha = 1) norm_post_NIG(x, prior_mu, prior_sigma, prior_alpha = prior_alpha)[3]

#########################
# Convenience functions #
#########################

lw = 1
lines <- function(...) geom_line(linewidth = 0.3 * lw, ...)
circles <- function(...) geom_point(pch = 21, fill = 'white', stroke = 0.5 * lw, ...)
ribbon <- function(...) geom_ribbon(col = NA, alpha = 0.2)

##############
### Ricing ###
##############

# Default theme (based on ggthemr)
source('~/Code/R_preamble/theme.R')

#' Adaptive palette (discrete)
#' from https://nanx.me/blog/post/ggplot2-color-interpolation/
#'
#' Create a discrete palette that will use the first `n` colors from
#' the supplied color values when the palette has enough colors.
#' Otherwise, use an interpolated color palette.
#'
#' @param values Color values.
pal_ramp <- function(values, use_raw = TRUE) {
  force(values)
  function(n) {
    if (use_raw) {
      values[seq_len(n)]
    } else {
      colorRampPalette(values, alpha = TRUE)(n)
    }
  }
}

#' Adaptive color palette generator for ggsci color palettes using `pal_ramp()`.
#'
#' @param name Color palette name in ggsci
#' @param palette Color palette type in ggsci
#' @param alpha Transparency level, a real number in (0, 1].
#'
#' @details See `names(ggsci:::ggsci_db)` for all color palette names in ggsci.
#' See `names(ggsci:::ggsci_db$"pal")` for available palette types under
#' the palette `pal`.
pal_adaptive <- function(palette, alpha = 1, use_raw = FALSE) {
  if (alpha > 1L | alpha <= 0L) stop("alpha must be in (0, 1]")

  #raw_cols <- ggsci:::ggsci_db[[name]][[palette]]
  raw_cols <- palette
  raw_cols_rgb <- col2rgb(raw_cols)
  alpha_cols <- rgb(
    raw_cols_rgb[1L, ], raw_cols_rgb[2L, ], raw_cols_rgb[3L, ],
    alpha = alpha * 255L, names = names(raw_cols),
    maxColorValue = 255L
  )

  pal_ramp(unname(alpha_cols), use_raw = use_raw)
}

#' Adaptive color scales
#'
#' @inheritParams pal_adaptive
#' @param ... additional parameters for [ggplot2::discrete_scale()].

navalue = 'grey60'
scale_color_adaptive <- function(palette, alpha = 1, use_raw = FALSE, ...) {
  ggplot2::discrete_scale("colour", "adaptive", pal_adaptive(palette, alpha, use_raw = use_raw), na.value = navalue, ...)
}

scale_fill_adaptive <- function(palette, alpha = 1, use_raw = FALSE, ...) {
  ggplot2::discrete_scale("fill", "adaptive", pal_adaptive(palette, alpha, use_raw = use_raw), na.value = navalue, ...)
}

# Backup palette (e.g. for datetime)
scale_color_def = function(...){
    scale_color_gradientn(colors = wes_palette('Zissou1'), ...)
}

scale_fill_def = function(...){
    scale_fill_gradientn(colors = wes_palette('Zissou1'), ...)
}

# Set default colors
library(wesanderson)
def_col_d = function(){return (scale_color_adaptive(wes_palette('Zissou1')))}
def_col_c = function(){return (scale_color_gradientn(colors = wes_palette("Zissou1")))}
def_fill_d = function(){return (scale_fill_adaptive(wes_palette('Zissou1')))}
def_fill_c = function(){return (scale_fill_gradientn(colors = wes_palette("Zissou1")))}

# Manual access
Zissou1 = c("#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00")
pal = c("#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00")

options(
  ggplot2.discrete.colour = def_col_d, ggplot2.continuous.colour = def_col_c,
  ggplot2.discrete.fill = def_fill_d, ggplot2.continuous.fill = def_fill_c
)
