\name{fitsad}
\alias{fitsad}
\alias{fitbs}
\alias{fitexp}
\alias{fitgamma}
\alias{fitgeom}
\alias{fitlnorm}
\alias{fitls}
\alias{fitmzsm}
\alias{fitnbinom}
\alias{fitpareto}
\alias{fitpoilog}
\alias{fitpower}
\alias{fitpowbend}
\alias{fitvolkov}
\alias{fitweibull}

\title{ML fitting of species abundance distributions}

\description{Fits probability distributions for
  abundances of species in a sample or assemblage by maximum likelihood.
}

\usage{
fitsad(x, sad = c("bs","exp", "gamma","geom","lnorm","ls","mzsm","nbinom","pareto",
       "poilog","power", "powbend", "volkov","weibull"), \dots)

fitbs(x, trunc, \dots)

fitexp(x, trunc = NULL, start.value, \dots)

fitgamma(x, trunc, start.value,  \dots)

fitgeom(x, trunc = 0, start.value, \dots)

fitlnorm(x, trunc, start.value,  \dots)

fitls(x, trunc, start.value, upper = length(x), \dots)

fitmzsm(x, trunc, start.value, upper = length(x), \dots)

fitnbinom(x, trunc = 0, start.value, \dots)

fitpareto(x, trunc, start.value,  upper = 20, \dots)

fitpoilog(x, trunc = 0, \dots)

fitpower(x, trunc, start.value, upper = 20, \dots)

fitpowbend(x, trunc, start.value, \dots)

fitvolkov(x, trunc, start.value, \dots)

fitweibull(x, trunc, start.value,  \dots)
}

\arguments{
  \item{x}{vector of (positive integer) quantiles. In the context of
    SADs, some abundance measurement (e.g., number of individuals,
    biomass) of species in a sample or ecological assemblage.
  }
  \item{sad}{
    character; root name of community sad distribution to be fitted.
    \kbd{"bs"} for broken-stick distribution,
    \kbd{"exp"} for exponential distribution,
    \kbd{"gamma"} for gamma distribution,
    \kbd{"geom"} for geometric distributions (not geometric series rad model,
    \code{\link{dgs}}),
    \kbd{"lnorm"} for lognormal,
    \kbd{"ls"} for Fisher's log-series,
    \kbd{"mzsm"} for Alonso & McKane's neutral metacommunity distribution,
    \kbd{"nbinom"} for negative binomial,
    \kbd{"pareto"} for Pareto distribution,
    \kbd{"poilog"} for Poisson-lognormal distribution,
    \kbd{"power"} for power-law distribution,
    \kbd{"powbend"} for Pueyo's power-bend distribution,
    \kbd{"volkov"} for Volkov's et al. neutral community distribution,
    \kbd{"weibull"} for Weibull distribution.
  }

  \item{trunc}{
    non-negative integer, \code{trunc > min(x)}; truncation point to fit a
    truncated distribution. For "poilog", "geom" and "nbinom" set
    trunc=NULL to avoid zero-truncation (see details).
  }

  \item{start.value}{
    named numeric vector; starting values of free parameters to be
    passed to \code{mle2}.
    Parameters should be named as in the corresponding density
    function,
    and in the same order.
  }
  
  \item{upper}{
    real positive; upper bound for Brent's one-parameter optimization
    method (default), for fits that use this method by default. See
    details and \code{optim}.
  }
  \item{\dots}{
    in \code{fitsad} further arguments to be passed to the specific fitting
    function (most used are \code{trunc}, \code{start.value})
    In the specific fitting functions further
    arguments to be passed to \code{mle2}.
  }
}

\value{
  An object of \code{fitsad-class} which inherits from \code{mle2-class}
  and thus has methods for handling
  results of maximum likelihood fits from \code{mle2}
  and also specific methods to handle SADs models
  (see \code{\link{fitsad-class}}).
}

\details{ \code{fitsad} is simply a wrapper that calls the specific
  functions to fit the distribution chosen with the argument
  \code{sad}. Users can interchangeably use \code{fitsad} or the
  individual functions detailed below (e.g. \code{fitsad(x, sad="geom",
  \dots)} is the same as \code{fitgeom(x, \dots)} and so on).
  
  The distributions are fitted by the maximum likelihood method using
  numerical optimization, with \code{mle2}.  The resulting object is of
  \code{fitsad-class} which can be handled with \code{mle2} methods for
  fitted models and has also some additional methods for SADs models
  (see \code{\link{fitsad-class}} and examples).
  
  Functions \code{fitexp}, \code{fitgamma}, \code{fitlnorm} and
  \code{fitweibull} fit the standard continuous distributions most used
  as SADs models.  As species with null abundances in the sample are in
  general unknown, the fit to continuous distributions can be improved
  by truncation to some value above zero (see example). Convergence
  problems can occur when fitting with truncation, and can be solved
  with sensible starting values.  \code{fitgamma} uses Chapman's (1956)
  fitting method to find starting values for the truncated gamma
  distribution, and \code{fitweibull} uses Rinne's (2009, p. 413) method
  (thanks to Mario Jose Marques-Azevedo).

  Functions \code{fitgeom}, \code{fitnbinom} fits geometric and negative
  binomial distributions which are two discrete standard distributions
  also used to fit SADs. Since species with zero individuals in the
  sample are in general unknown, these functions fit by default
  zero-truncated distributions.  To avoid zero-truncation set
  \code{trunc=NULL}.  Using the geometric distribution as a SAD model is
  not to be confounded for fitting the Geometric series
  \code{\link{fitgs}} as a rank-abundance distribution (RAD) model.

  Function \code{fitls} implements the original numerical recipe by
  Fisher (1943) to fit the log-series distribution, given a vector of
  species abundances.  Alonso et al. (2008, supplementary material)
  showed that this recipe gives the maximum likelihood estimate of
  Fisher's alpha, the single parameter of the log-series.  Fitting is
  done through numerical optimization with the \code{uniroot} function,
  following the code of the function \code{fishers.alpha} of the
  \pkg{untb} package. After that, the estimated value of alpha parameter
  is used as the starting value to get the Log-likelihood from the
  log-series density function \code{dls}, using the function
  \code{mle2}.  The total number of individuals in the sample, \code{N},
  is treated as a fixed parameter in this implementation, in order to
  maintain coherence with similar parameters from \code{fitvolkov} and
  \code{fitmzsm} (see below).  Fixed parameters in the model
  specification do not contribute to the model degrees of freedom, and
  are not accounted in standard error calculations.

  Functions \code{fitpower}, \code{fitpowbend} and \code{fitpareto} fit
  power-law distributions with one and two-parameters, that have been
  suggested as SADs models.  The implementations of power and power-bend
  are for discrete distributions that do not include zeroes. The Pareto
  distribution is continuous and includes all non-negative numbers.
  Fisher's logseries are a special case of the power-bend, see
  \code{\link{dpowbend}} and Pueyo (2006).
  
  Function \code{fitbs} fits the Broken-stick distribution (MacArthur
  1960). It is defined only by the observed number of elements \code{S}
  in the collection and collection size \code{N}.  Thus once a sample is
  taken, the Broken-stick has no free parameters.  Therefore, there is
  no actual fitting, but still \code{fitbs} calls \code{mle2} with fixed
  parameters \code{N} and \code{S} and \code{eval.only=TRUE} to return
  an object of class \code{fitsad} to keep compatibility with other SADs
  models fitted to the same data.  Therefore, the resulting objects
  allows most of the operations with SAD models, such as comparison with
  other models through model selection, diagnostic plots and so on (see
  \code{\link{fitsad-class}}).

  Function \code{fitpoilog} fits the Poisson-lognormal distribution.
  This is a compound distributions that describes the abundances of
  species in Poisson sample of community that follows a lognormal
  SAD. This is a sampling model of SAD, which is approximated by the
  \sQuote{veil line} truncation of the lognormal (Preston 1948). The
  Poisson-lognormal is the analytic solution for this sampling model, as
  Fisher's log-series is a analytic limit case for a Poisson-gamma
  (a.k.a negative binomial) distribution. As geometric and negative
  binomial distributions, the Poisson-lognormal includes zero, but the
  fit is zero-truncated by default, as for \code{fitgeom},
  \code{fitnbinom}.  To avoid zero-truncation set \code{trunc=NULL}.

  \code{fitmzsm} fits the metacommunity Zero-sum multinomial
  distribution \code{\link{dmzsm}} from the Neutral Theory of
  Biodiversity (Alonso and McKane 2004). The mZSM describes the SAD of a
  sample taken from a neutral metacommunity under random drift.  It has
  two parameters, the number of individuals in the sample \code{J} and
  \code{theta}, the \sQuote{fundamental biodiversity number}.  Because
  \code{J} is known from the sample size, the fit resumes to estimate a
  single parameter, \code{theta}. By default, \code{fitmzsm} fits mZSM
  to a vector of abundances with Brent's one-dimensional method of
  optimization (see \code{optim}). The log-series distribution (Fisher
  et al. 1943) is a limiting case of mZSM (Hubbel 2001), and
  \code{theta} tends to Fisher's alpha as \code{J} increases. In
  practice the two models provide very similar fits to SADs (see
  example).

  Function \code{fitvolkov} fits the SAD model for a community under
  neutral drift with immigration (Volkov et al. 2003).  The model is a
  stationary distribution deduced from a stochastic process compatible
  with the Neutral Theory of Biodiversity (Hubbell 2001). It has two
  free parameters, the \sQuote{fundamental biodiversity number}
  \code{theta}, and the immigration rate \code{m} (see
  \code{\link{dvolkov}}) \code{fitvolkov} builds on function
  \code{volkov} from package \pkg{untb} to fit Volkov's \emph{et al.}
  SAD model to a vector of abundances.  The fit can be extremely slow
  even for vectors of moderate size.  }

\references{

  Alonso, D. and McKane, A.J.  2004.  Sampling Hubbell's neutral
  model of biodiversity. \emph{Ecology Letters 7}:901--910

  Alonso, D. and Ostling, A., and Etienne, R.S. 2008 The implicit
  assumption of symmetry and the species abundance
  distribution. \emph{Ecology Letters, 11}: 93-105.

  Bolker, B. and R Development Core Team 2012. bbmle: Tools for
  general maximum likelihood estimation. R package version 1.0.5.2.
  http://CRAN.R-project.org/package=bbmle

  Chapman, D. G. 1956. Estimating the parameters of a truncated gamma
  distribution. \emph{The Annals of Mathematical Statistics, 27(2)}:
  498--506.

  Fisher, R.A, Corbert, A.S. and Williams, C.B. (1943) The Relation
  between the number of species and the number of individuals in a
  random sample of an animal population.
  \emph{The Journal of Animal Ecology, 12(1)}: 42--58.

  \enc{Grøtan}{Grotan}, V. and Engen, S. 2008. poilog: Poisson lognormal and
  bivariate Poisson lognormal distribution. R package version 0.4.
  
  Hankin, R.K.S. 2007. Introducing untb, an \R Package For Simulating Ecological Drift
  Under the Unified Neutral Theory of Biodiversity. \emph{Journal
    of Statistical Software 22} (12). 
  
  Hubbell, S.P.  2001.  \emph{The Unified Neutral Theory of Biodiversity}.
  Princeton University Press
  
  MacArthur, R.H. 1960. On the relative abundance of species.
  \emph{Am Nat 94}:25--36.
  
  Magurran, A.E. 1989. \emph{Ecological diversity and its measurement}.
  Princenton University Press.
  
  Preston, F.W. 1948. The commonness and rarity of species.
  \emph{Ecology 29}: 254--283.

  Pueyo, S. (2006) \emph{Diversity: Between neutrality and structure},
  Oikos 112: 392-405.

  Rinne, H. 2009. \emph{The Weibull distribution: a handbook}. CRC Press

  Volkov, I.,  Banavar, J. R., Hubbell, S. P., Maritan, A. 2003.
  Neutral theory and relative species abundance in ecology.  \emph{Nature 424}:1035--1037.
}

\source{all fitting functions builds on \code{mle2} and methods
  from \pkg{bbmle} package (Bolker 2012), which in turn builds on
  \code{mle} function and associated classes and methods;
  \code{fitls} and \code{fitvolkov} use codes and functions from
  \pkg{untb} package (Hankin 2007); \code{fitpoilog} builds on
  \pkg{poilog} package (\enc{Grøtan}{Grotan} & Engen 2008).
}

\seealso{ \code{\link{dls}}, \code{\link{dmzsm}}, \code{\link{dpareto}},
  \code{\link{dpoilog}}, \code{\link{dpower}},\code{ \link{dvolkov}},
  \code{\link{dpowbend}} for corresponding density functions created for
  fitting SADs; standard distributions \code{dexp}, \code{dgamma},
  \code{dgeom}, \code{dlnorm}, \code{dnbinom}, \code{dweibull};
  \code{\link{fitsad-class}}.  }

\examples{

## Magurran (1989) example 5:
## birds in an Australian forest
mag5 <- c(103, 115, 13, 2, 67, 36, 51, 8, 6, 61, 10, 21,
          7, 65, 4, 49, 92, 37, 16, 6, 23, 9, 2, 6, 5, 4, 
          1, 3, 1, 9, 2)
mag5.bs <- fitsad(mag5, "bs") 
summary(mag5.bs)## no estimated coefficient
coef(mag5.bs) ## fixed coefficients N and S
## Diagnostic plots
par(mfrow=c(2, 2))
plot(mag5.bs)
par(mfrow=c(1, 1))


data(moths) #Fisher's moths data
moths.mzsm <- fitmzsm(moths) ## same as fitsad(moths, sad="mzsm")
## fit to log-series
moths.ls <- fitsad(moths, sad="ls")
coef(moths.ls)
coef(moths.mzsm) ## Compare with theta=38.9, Alonso&McKanne (2004)
## Diagnostic plots
par(mfrow=c(2, 2))
plot(moths.mzsm)
par(mfrow=c(1, 1))
## Graphical comparison
plot(rad(moths))
lines(radpred(moths.ls))
lines(radpred(moths.mzsm), col="red", lty=2)
legend("topright", c("log-series", "mZSM"), lty=1, col=c("blue","red"))
## Two more models: truncated lognormal and Poisson-lognormal
moths.ln <- fitsad(moths, "lnorm", trunc=0.5)
moths.pln <- fitsad(moths, "poilog")
## Model selection
AICtab(moths.ln, moths.pln, moths.ls, moths.mzsm, weights=TRUE)

## Biomass as abundance variable
data(ARN82.eB.apr77) #benthonic marine animals
AR.ln <- fitsad(ARN82.eB.apr77, sad="lnorm")
AR.g <- fitsad(ARN82.eB.apr77, sad="gamma")
AR.wb <- fitsad(ARN82.eB.apr77, sad="weibull")
plot(octav(ARN82.eB.apr77))
lines(octavpred(AR.ln))
lines(octavpred(AR.g), col="red")
lines(octavpred(AR.wb), col="green")
legend("topright", c("lognormal", "gamma", "weibull"),lty=1, col=c("blue","red", "green"))
AICctab(AR.ln, AR.g, AR.wb, nobs=length(ARN82.eB.apr77), weights=TRUE) 
}

\author{Paulo I Prado \email{prado@ib.usp.br}, Murilo Dantas Miranda and
  Andre Chalom,
  after Ben Bolker, R Core Team, Robin Hanking, Vidar \enc{Grøtan}{Grotan} and Steinar Engen.}

\encoding{UTF-8}
