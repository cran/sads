\name{dbs}

\alias{dbs}
\alias{pbs}
\alias{qbs}
\alias{drbs}
\alias{prbs}
\alias{qrbs}
\alias{rrbs}
\alias{rbs}

\title{MacArthur's Broken-stick distribution}

\description{
Density, distribution function, quantile function and random generation for
the Broken-stick distribution with parameters \code{N} and \code{S}.
}

\usage{
dbs( x, N, S, log = FALSE )
pbs( q, N, S, lower.tail = TRUE, log.p = FALSE )
qbs( p, N, S, lower.tail = TRUE, log.p = FALSE )
rbs( n, N, S )
drbs( x, N, S, log = FALSE )
prbs( q, N, S, lower.tail = TRUE, log.p = FALSE )
qrbs( p, N, S, lower.tail = TRUE, log.p = FALSE )
rrbs( n, N, S)
}

\arguments{
  \item{x}{
    vector of (non-negative integer) quantiles. In the context of
    species abundance distributions, this is a vector of abundances (for
    \code{dbs}) or abundance ranks (for \code{drbs})
    of species in a sample.
  }
  \item{q}{vector of (non-negative integer) quantiles. In the context of
    species abundance distributions, a vector of abundances
    (for \code{dbs})
    or abundance ranks (for \code{drbs}) of species in a sample.
  }
  \item{n}{number of random values to return.}
  \item{p}{
    vector of probabilities.
  }
  \item{N}{
    positive integer 0 < N < Inf, sample size. In the context of
    species abundance distributions, the sum of abundances of individuals in a sample.
  }
  \item{S}{
    positive integer 0 < S < Inf, number of elements in a collection.
    In the context of species abundance distributions, the number of species in a sample.
  }
  \item{log, log.p}{
    logical; if TRUE, probabilities p are given as log(p).
  }
  \item{lower.tail}{
    logical; if TRUE (default), probabilities are P[X <= x],
    otherwise, P[X > x].
  }
}

\details{
  
  The Broken-stick distribution was proposed as a model for
  the expected abundance of elements in a collection:

  \deqn{n(i) = \frac{N}{S} \sum_{k=i}^S 1/k}{n(i) = N/S (sum(from k=i to S) 1/k)}

  where n(i) is the abundance in the i-th most abundant element
  (MacArthur 1960, May 1975).
  Hence the probability (or expected proportion of occurrences)
  in the i-th element is

  \deqn{p(i) = \frac{n(i)}{S} = S^{-1}\sum_{k=i}^S 1/k}{p(i) = n(i)/N = (sum(from k=i to S) 1/k) / S}

  \code{[dpq]rbs} stands for "rank-abundance Broken-stick" and return
  probabilities and quantiles based on the expression above, for p(i).
  Therefore, \code{[dpq]rbs} can be used as a rank-abundance model
  for species' ranks in a sample or in a biological community
  see \code{\link{fitrad}}.

  The probability density for a given abundance value in the Broken-stick model
  is given by 

  \deqn{p(x) = \frac{S-1}{N} \left( 1 -  \frac{x}{N} \right)^{S-2}}{p(x) = (1 - x/N)^(S-2) (S - 1)/N}

  Where x is the abundance of a given element in the collection (May 1975).
  \code{[dpq]bs} return probabilities and quantiles according to the
  expression above for p(x).
  Therefore, \code{[dpq]bs} can be used as a
  species abundance model
  see \code{\link{fitsad}}.
  }

\value{
  \code{dbs}  gives the (log) density and  \code{pbs} gives the (log)
  distribution function of abundances,  and \code{qbs} gives the
  corresponding quantile function.
  \code{drbs}  gives the (log) density and  \code{prbs} gives the (log)
  distribution function of ranks,  and \code{qrbs} gives the
  corresponding quantile function.
}

\references{

  MacArthur, R.H. 1960. On the relative abundance of species.
  \emph{Am Nat 94}:25--36.
  
  May, R.M. 1975. Patterns of Species Abundance and Diversity. In
  Cody, M.L. and Diamond, J.M. (Eds) \emph{Ecology and Evolution of
  Communities}. Harvard University Press. pp 81--120.
  }

\author{Paulo I Prado \email{prado@ib.usp.br} and Murilo Dantas Miranda.}

\seealso{
  \code{\link{fitbs}} and \code{\link{fitrbs}} to fit the Broken-stick distribution
  as a abundance (SAD) and rank-abundance (RAD) model.
}

\examples{
x <- 1:25
PDF <- drbs(x=x, N=100, S=25)
CDF <- prbs(q=x, N=100, S=25)
par(mfrow=c(1,2))
plot(x,CDF, ylab="Cumulative Probability", type="b",
     main="Broken-stick rank distribution, CDF")
plot(x,PDF, ylab="Probability", type="h",
     main="Broken-stick rank distribution, PDF")
par(mfrow=c(1,1))

## quantile is the inverse of CDF
all.equal( qrbs( CDF, N=100, S=25), x) # should be TRUE
}
