# dispfit

R package to estimate species dispersal kernels

The *dispfit* package is an R software application developed to intuitively and comprehensively estimate dispersal kernels from dispersal data. *dispfit* fits and compares different families of parameterized functions to describe and predict dispersal distances. It includes 9 well-known and commonly used distributions, computing goodness-of-fit and model selection statistics, and estimating each distribution's parameters, along with their first four moments (mean, standard deviation, skewness, and kurtosis).

To install *dispfit,* run the following commands in R:

    # Install 'devtools' package, if needed

    install.packages("devtools") 

    devtools::install_github("https://github.com/apferreira/dispfit")

    library(dispfit)
