
# context("RABR with continuous endpoints")
library(RABR)
library(parallel)
library(doParallel)

test_that("Length of unadjusted probabilities is the length of number of active treatment groups", {
  RABRtest = RABRcontinuous(
    MeanVec = c(0.15, 0.28, 2),
    SdVec = c(1, 1, 1),
    M = 100,
    N = 200,
    R = c(2, 1, 1),
    Nitt = 10,
    Alpha = 0.025,
    Ncluster = 1,
    Seed = 1,
    MultiMethod = "bonferroni")

  expect_equal(length(RABRtest$ProbUnadj), 2)
})

test_that("Length of adjusted probabilities is the length of number of active treatment groups", {
  RABRtest = RABRcontinuous(
    MeanVec = c(0.15, 0.28, 2),
    SdVec = c(1, 1, 1),
    M = 100,
    N = 200,
    R = c(2, 1, 1),
    Nitt = 10,
    Alpha = 0.025,
    Ncluster = 1,
    Seed = 1,
    MultiMethod = "bonferroni")

  expect_equal(length(RABRtest$ProbAdj), 2)
})

test_that("Length of selected probabilities is the length of number of active treatment groups", {
  RABRtest = RABRcontinuous(
    MeanVec = c(0.15, 0.28, 2),
    SdVec = c(1, 1, 1),
    M = 100,
    N = 200,
    R = c(2, 1, 1),
    Nitt = 10,
    Alpha = 0.025,
    Ncluster = 1,
    Seed = 1,
    MultiMethod = "bonferroni")

  expect_equal(length(RABRtest$ProbAdjSelected), 2)
})

test_that("Length of overall probability is 1", {
  RABRtest = RABRcontinuous(
    MeanVec = c(0.15, 0.28, 2),
    SdVec = c(1, 1, 1),
    M = 100,
    N = 200,
    R = c(2, 1, 1),
    Nitt = 10,
    Alpha = 0.025,
    Ncluster = 1,
    Seed = 1,
    MultiMethod = "bonferroni")

  expect_equal(length(RABRtest$ProbAdjOverall), 1)
})

test_that("Length of ASN is the number of treatment groups", {
  RABRtest = RABRcontinuous(
    MeanVec = c(0.15, 0.28, 2),
    SdVec = c(1, 1, 1),
    M = 100,
    N = 200,
    R = c(2, 1, 1),
    Nitt = 10,
    Alpha = 0.025,
    Ncluster = 1,
    Seed = 1,
    MultiMethod = "bonferroni")

  expect_equal(length(RABRtest$ASN), 3)
})














