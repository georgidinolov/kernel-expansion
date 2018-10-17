#!/bin/bash

bazel --output_user_root=/home/gdinolov-tmp build //src/brownian-motion:generate-likelihood-points-with-mle

bazel --output_user_root=/home/gdinolov-tmp build //src/finite-element-igraph:kriging-test

./bazel-bin/src/brownian-motion/generate-likelihood-points-with-mle 10 80 ./src/kernel-expansion/documentation/chapter-4/spy-ftse.csv ./src/kernel-expansion/documentation/chapter-4/spy-ftse-interpolator.csv

./bazel-bin/src/finite-element-igraph/kriging-test 50 0.00 0.40 0.40 0.001 ./src/kernel-expansion/documentation/chapter-4/spy-ftse-interpolator-dx-1200 ./src/kernel-expansion/documentation/chapter-4/spy-ftse-interpolator.csv


# squares <- function(x, y, tt) {
#     x[6] = 0
#     out = sum((y-exp(-tt*x[1])*(x[2]*tt^4 + x[3]*tt^3 + x[4]*tt^2 + x[5]*tt + x[6]))^2)
#     return (out)
# }



# x
# [1]   0.26468530  -0.03333049  -0.44839950   4.18190921 -31.55219702
# [6]  47.59188573