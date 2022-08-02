
#' Simulate RABR for continuous endpoints to evaluate operating characteristics
#'
#' @importFrom stats t.test
#' @importFrom stats aov
#' @importFrom stats rmultinom
#' @importFrom stats rnorm
#' @importFrom stats runif
#' @importFrom stats sd
#' @importFrom foreach %dopar%
#' @importFrom stats p.adjust
#' @importFrom stats prop.test
#'
#' @param MeanVec Vector of response mean for placebo and active treatment groups.
#' @param SdVec Vector of standard deviation for placebo and active treatment groups.
#' @param M Total sample size of burn-in period.
#' @param N Total sample size of RABR. Must be larger than M.
#' @param R Randomization vector for placebo and active treatment groups.
#' @param Nitt Number of simulation iterations.
#' @param Alpha One-sided significance level.
#' @param Ncluster Number of clusters for parallel computing.
#' @param Seed Random seed.
#' @param MultiMethod Multiplicity adjustment method. Must be one of the following values "holm", "hochberg", "hommel", "bonferroni", or "dunnett".
#'
#' @details The \code{MeanVec} is a vector of response mean for placebo and active treatment groups, while \code{SdVec} is for standard deviation. They should be with the same length. The current package supports 2 or 3 active treatment groups. Note that a larger response corresponds to a better outcome.
#' @details The \code{M} is the total sample size of burn-in period with equal randomization. The total sample size \code{N} should be larger than N. The choice of \code{M} can be selected by comparing simulations from several candidate values. The \code{R} is a pre-specified randomization vector, where the first element is for placebo, and the next one for the best performing group, up to the worst performing group.
#' @details The \code{Alpha} is the one-sided significance level. The \code{MultiMethod} can be set at "holm" for Holm, "hochberg" for Hochberg, "hommel" for Hommel, "bonferroni" for Bonferroni, or "dunnett" for Dunnett procedures.
#'
#' @return ProbUnadj: Probability of rejecting each elementary null hypothesis without multiplicity adjustment
#' @return ProbAdj: Probability of rejecting each elementary null hypothesis with multiplicity adjustment
#' @return ProbAdjSelected: Probability of selecting and confirming the efficacy of each active treatment group
#' @return ProbAdjOverall: Probability of rejecting at least one elementary null hypothesis with multiplicity adjustment
#' @return ASN: Average sample size of placebo and active treatment groups
#' @export
#' @author Tianyu Zhan (tianyu.zhan.stats@gmail.com)
#' @references Zhan, T., Cui, L., Geng, Z., Zhang, L., Gu, Y., & Chan, I. S. (2021). A practical response adaptive block randomization (RABR) design with analytic type I error protection. Statistics in Medicine, 40(23), 4947-4960.
#' @references Cui, L., Zhan, T., Zhang, L., Geng, Z., Gu, Y., & Chan, I. S. (2021). An automation-based adaptive seamless design for dose selection and confirmation with improved power and efficiency. Statistical Methods in Medical Research, 30(4), 1013-1025.
#' @examples ## Consider an example with three active treatment
#' @examples ## groups and a placebo. Suppose that the response
#' @examples ## mean for placebo is 0.43 and 0.48, 0.63, and 1.2
#' @examples ## for three active treatment groups. The standard
#' @examples ## deviation is 1 for all groups. The total sample
#' @examples ## size is N = 120 with a burn-in period M = 60. We
#' @examples ## use the randomization vector of (8, 9, 2, 1),
#' @examples ## which means that placebo, the best performing
#' @examples ## group, the second-best group, and the worst group
#' @examples ## have randomization probabilities 8/20, 9/20, 2/20
#' @examples ## 1/20, respectively. The one-sided significance
#' @examples ## level is considered at 2.5%. Nitt = 100 is for
#' @examples ## demonstration, and should be increased to 10^5
#' @examples ## in practice.
#' @examples ##
#' @examples library(parallel)
#' @examples library(doParallel)
#' @examples RABR.fit = RABRcontinuous(
#' @examples            MeanVec = c(0.43, 0.48, 0.63, 1.2),
#' @examples            SdVec = c(1, 1, 1, 1),
#' @examples            M = 60,
#' @examples            N = 120,
#' @examples            R = c(8, 9, 2, 1),
#' @examples            Nitt = 100,
#' @examples            Alpha = 0.025,
#' @examples            Ncluster = 2,
#' @examples            Seed = 12345,
#' @examples            MultiMethod = "dunnett")
#' @examples ##
#' @examples ## Probability of rejecting each elementary null
#' @examples ## hypothesis without multiplicity adjustment
#' @examples    print(RABR.fit$ProbUnadj)
#' @examples ##
#' @examples ## Probability of rejecting each elementary null
#' @examples ## hypothesis with multiplicity adjustment
#' @examples    print(RABR.fit$ProbAdj)
#' @examples ##
#' @examples ## Probability of selecting and confirming the
#' @examples ## efficacy of each active treatment group
#' @examples    print(RABR.fit$ProbAdjSelected)
#' @examples ##
#' @examples ## ProbAdjOverall Probability of rejecting at
#' @examples ## least one elementary null hypothesis
#' @examples ## with multiplicity adjustment
#' @examples    print(RABR.fit$ProbAdjOverall)
#' @examples ##
#' @examples ## ASN Average sample size of placebo and active
#' @examples ## treatment groups
#' @examples    print(RABR.fit$ASN)
#'
#'
#'
#'
RABRcontinuous = function(MeanVec, SdVec, M, N, R, Nitt, Alpha, Ncluster = 1, Seed = 12345, MultiMethod){

  ## check input
  check.input.ind.vec = rep(0, 13)
  check.input.message = "Please check the following error message:"
  if (!length(MeanVec)==length(SdVec)|!length(MeanVec)==length(R)){
    check.input.ind.vec[1] = 1
    check.input.message = paste0(check.input.message, "\n",
                                 "The lengths of MeanVec, SdVec and R should be the same.")
  }
  if (length(MeanVec)>4|length(MeanVec)<3){
    check.input.ind.vec[2] = 1
    check.input.message = paste0(check.input.message, "\n",
                                 "The number of active treatment groups should be either 2 or 3.")
  }
  if (sum(SdVec<=0)>0){
    check.input.ind.vec[3] = 1
    check.input.message = paste0(check.input.message, "\n",
                                 "The standard deviation should be positive.")
  }
  if (M>=N){
    check.input.ind.vec[4] = 1
    check.input.message = paste0(check.input.message, "\n",
                                 "M should be smaller than N.")
  }
  if (M<2*length(MeanVec)){
    check.input.ind.vec[5] = 1
    check.input.message = paste0(check.input.message, "\n",
                                 "M should not be smaller than 2 * number of treatment groups.")
  }
  if (sum(R<0)>0){
    check.input.ind.vec[6] = 1
    check.input.message = paste0(check.input.message, "\n",
                                 "R should be non-negative.")
  }
  if (sum(R)==0){
    check.input.ind.vec[7] = 1
    check.input.message = paste0(check.input.message, "\n",
                                 "At least one element in R should be positive.")
  }
  if (sum(diff(R)[-1]>0)>0){
    check.input.ind.vec[8] = 1
    check.input.message = paste0(check.input.message, "\n",
                                 "Elements in R exclusing the first one for placebo should be in a non-increasing order.")
  }
  if (Nitt<=0 | !Nitt==round(Nitt)){
    check.input.ind.vec[9] = 1
    check.input.message = paste0(check.input.message, "\n",
                                 "Nitt should be a positive integer.")
  }
  if (Alpha>=1 | Alpha<=0){
    check.input.ind.vec[10] = 1
    check.input.message = paste0(check.input.message, "\n",
                                 "Alpha should be between 0 and 1.")
  }
  if (Ncluster<=0 | !Ncluster==round(Ncluster)){
    check.input.ind.vec[11] = 1
    check.input.message = paste0(check.input.message, "\n",
                                 "Ncluster should be a positive integer.")
  }
  if (!is.numeric(Seed)){
    check.input.ind.vec[12] = 1
    check.input.message = paste0(check.input.message, "\n",
                                 "Seed should be numeric.")
  }
  if ((sum(MultiMethod==c("holm", "hochberg", "hommel", "bonferroni", "dunnett"))==0)){
    check.input.ind.vec[13] = 1
    check.input.message = paste0(check.input.message, "\n",
                                 "MultiMethod should be one the following methods: holm, hochberg, hommel, bonferroni, or dunnett.")
  }

  if (sum(check.input.ind.vec)>0){
    stop(check.input.message)
  }


  resp.vec = MeanVec
  sd.vec = SdVec
  m.total = M
  n.arm = length(resp.vec)
  n.active.arm = n.arm -1
  n.total = N
  n.look = n.total-m.total
  block = R
  n.itt = Nitt
  one.sided.alpha = Alpha
  n.cluster = Ncluster
  seed.in = Seed

RAR.dunnett.func = function(current.sample.func, n.arm.func){

  data.dunnett.stage.1 = data.frame("x" = unlist(current.sample.func),
                                    "group" = as.factor((
                                      unlist(sapply(1:n.arm.func, function(x){rep(x, length(unlist(current.sample.func[[x]])))}))
                                    )))

  dunnett.stage.1.anova = aov(x ~ group, data.dunnett.stage.1)
  dunnett.stage.1.fit = multcomp::glht(dunnett.stage.1.anova,
                                       linfct = multcomp::mcp(group = "Dunnett"),
                             alternative = "greater")
  dunnett.stage.1.summary = summary(dunnett.stage.1.fit,
                                    test=(multcomp::adjusted(type = "free")))

  p.value.dunnett.vec = dunnett.stage.1.summary$test$pvalues
  p.value.dunnett.vec = p.value.dunnett.vec + runif(n.arm.func-1, 0, 1)*10^(-6)
  return(p.value.dunnett.vec)
}

RAR.t.func = function(current.sample.list, n.active.arm){
  pbo.vec = current.sample.list[[1]]
  dec.vec.out = rep(NA, n.active.arm)
  for (i in 1:n.active.arm){
    trt.vec = current.sample.list[[1+i]]

    p.out = t.test(trt.vec, pbo.vec, alternative = "greater")$p.value

    dec.vec.out[i] = p.out
  }
  return(dec.vec.out)
}

cl <- parallel::makeCluster(n.cluster)
doParallel::registerDoParallel(cl)
pred = foreach::foreach(itt = 1:n.itt, .packages=c("multcomp")) %dopar% {

  set.seed(seed.in + itt)

  rand.initial.m.vec = rep(1/n.arm, n.arm)
  n.initial.m.vec = round(rand.initial.m.vec*m.total)

  current.sample.list = lapply(1:n.arm,
          function(x){rnorm(n.initial.m.vec[x], resp.vec[x], sd.vec[x])})

  for (i in 1:n.look){

    stand.mean.vec = sapply(2:n.arm, function(x){(mean(current.sample.list[[x]]))*
        sqrt(length(current.sample.list[[x]]))/sd(current.sample.list[[x]])})

    ## break tie
    if (!(length(unique(stand.mean.vec))==n.active.arm)){
      stand.mean.vec = stand.mean.vec + (n.active.arm:1)*0.00001
    }

    RAR.ratio.unadj = c(block[1],
              block[-1][match(stand.mean.vec,sort(stand.mean.vec, decreasing = TRUE))])

    new.sample.n = as.vector(rmultinom(1, 1, RAR.ratio.unadj/sum(RAR.ratio.unadj)))

    new.sample.list = lapply(1:n.arm,
                             function(x){rnorm(new.sample.n[x], resp.vec[x], sd.vec[x])})

    current.sample.list = mapply(c, current.sample.list, new.sample.list, SIMPLIFY=FALSE)

  }

  p.value.unadj.vec = RAR.t.func(current.sample.list, n.active.arm)

  if (MultiMethod=="dunnett"){
    dunnett.fit = RAR.dunnett.func(current.sample.list, n.arm)
    p.value.adj.vec = dunnett.fit
  } else {
    p.value.adj.vec = stats::p.adjust(p.value.unadj.vec, method = MultiMethod)
  }

  RAR.selected.arm = which.min(p.value.unadj.vec)

  dec.vec = rep(0, n.active.arm)
  dec.vec[RAR.selected.arm] = as.numeric(p.value.adj.vec[RAR.selected.arm]<=one.sided.alpha)

  ## ASN
  ASN.vec = sapply(current.sample.list, length)
  ASN.trt.vec = ASN.vec[-1]
  ASN.trt.s.vec = ASN.trt.vec[order(p.value.adj.vec, decreasing = FALSE)]

  clustlist = list("p.value.unadj.vec" = p.value.unadj.vec,
                   "p.value.adj.vec" = p.value.adj.vec,
                   "dec.vec" = dec.vec,
                   "ASN" = c(ASN.vec[1], ASN.trt.s.vec)
  )
  return(clustlist)


}

parallel::stopCluster(cl)

p.unadj.mat = p.adj.mat = dec.mat = matrix(NA, nrow = n.itt, ncol = n.active.arm)
ASN.mat = matrix(NA, nrow = n.itt, ncol = n.arm)

for (itt in 1:n.itt){
     pred.temp = pred[[itt]]
     p.unadj.mat[itt, ] = pred.temp$p.value.unadj.vec
     p.adj.mat[itt, ] = pred.temp$p.value.adj.vec
     ASN.mat[itt, ] = pred.temp$ASN
     dec.mat[itt, ] = pred.temp$dec.vec
}

power.unadj.output = apply(p.unadj.mat, 2, function(x){mean(x<=one.sided.alpha)})
power.adj.output = apply(p.adj.mat, 2, function(x){mean(x<=one.sided.alpha)})
power.adj.overal = 1-mean(apply(p.adj.mat, 1, function(x){sum(x>one.sided.alpha)})==n.active.arm)
ASN.output = apply(ASN.mat, 2, mean)
dec.output = apply(dec.mat, 2, mean)

return(list("ProbUnadj" = power.unadj.output,
            "ProbAdj" = power.adj.output,
            "ProbAdjSelected" = dec.output,
            "ProbAdjOverall" = power.adj.overal,
            "ASN" = ASN.output))

}























