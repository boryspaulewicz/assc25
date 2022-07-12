## Once you have the remotes package installed you can install the bhsdtr2 package like so:
##
## library(remotes)
## remotes::install_github('boryspaulewicz/bhsdtr2')

library(bhsdtr2)
library(ggplot2)

####################################################################################################
## SDT as an extremely simplified decision process model
##
## The simulation below explicitly represents the causal process assumed in SDT

n = 1000
## This is a series of randomly assigned stimuli
stimulus = rbinom(n, 1, .5)
## We will encode them as -1 or 1 instead of 0 and 1
stimulus = 2 * (stimulus - .5)
## because this way we can simplify the simulation
sensitivity = 1.5
## Notice that internal evidence samples depend on the stimulus class ...
evidence = rnorm(n, sensitivity * .5 * stimulus)
## ... but the decision criterion does not
criterion = 0.8
## This is how the motor response is generated according to the most common version of SDT
response = as.numeric(evidence > criterion)

## There is just one "participant" here so we can use a non-hierarchical model to estimate the SDT
## parameters. In this case, the easiest way is to simply fit a probit regression model like so:
(m1 = glm(response ~ I(stimulus / 2), family = binomial(link = probit)))
## Notice that (Intercept) = -criterion and the slope is sensitivity

## You can do the same thing with the bhsdtr2 package but it will take much more time
m2 = bhsdtr(c(dprim ~ 1, thr ~ 1), response ~ stimulus, data.frame(response, stimulus))
samples(m2, 'dprim')
samples(m2, 'thr')

####################################################################################################
## SDT for rating experiments

## We did not plot the model fit nor did we evaluate it in any other way. That is because the model
## is saturated and so it has to fit perfectly. However, we never know a priori how bad an SDT model
## is in describing our data so when possible it is a good idea to collect e.g., confidence ratings
## and see if an SDT model fits. The sim.sdt function in the bhsdtr2 package generates samples from
## SDT models with one or more threshold (~ decision criterion) and with possibly unequal
## variance. We will use the default values: d' = 1.5, thresholds = -2.1, -1.4, -0.7, 0, 0.7, 1.4,
## 2.1, variance of both evidence distributions is the same and fixed at 1
d = sim.sdt(n)
## We cannot fit this model using the glm function as before
m3 = bhsdtr(c(dprim ~ 1, thr ~ 1), r ~ stim, d)
samples(m3, 'dprim')
samples(m3, 'thr')
## It makes sense to see how the model fits
plot(m3)
## Notice that there are no error bars in this plot. That is because without the method='stan'
## argument the bhsdtr2 function only gives point estimates (joint maximum aposteriori point
## estimates to be exact) which is much faster than sampling from the posterior.

## Let's do some statistical inference
m4 = bhsdtr(c(dprim ~ 1, thr ~ 1), r ~ stim, d, method = 'stan')
## Now we can see the error bars but it does take time to even generate the plot
plot(m4)
## Unsurprisingly, there seems to be no evidence to reject the model
samples(m4, 'dprim')
samples(m4, 'thr')

## When printed the object returned by the samples function only shows a simple summary consisting
## of posterior means and the total number of posterior samples (the chains are glued toghether to
## make one chain)
s.dprim = samples(m4, 'dprim')
s.thr = samples(m4, 'thr')
## The bhsdtr_samples object is a three-dimensional array.
dim(s.dprim)
## The first dimension is the sample number, the second dimension is the dimension of the parameter
## - there is just one d' here so this dimension equals 1, and the third dimension is the condition,
## if any - there are no separate conditions here which means that there is just one condition.
dimnames(s.dprim)
## Notice that the second dimension is not equal to 1 for the thresholds
dim(s.thr)

## Let's fit a hierarchical version of SDT to the data that come with the bhsdtr2 package. This is
## how this dataset looks:
head(gabor)
## It is supposed to be an example of a typical dataset that an SDT model could be used for: there
## are repeated measurement, two classes of stimuli, and there are two regular predictors (order,
## which is a between-person factor, and duration, which is within=person) that may be associated
## with d' or the thresholds.

## Usually, you will not have the combined response in your data. To obtain one you can use the
## combined.response function
gabor$r = combined.response(gabor$stim, gabor$rating, gabor$acc)
## Let's do a quick fit first
m6 = bhsdtr(c(dprim ~ duration * order + (duration | id), thr ~ order + (1 | id)),
            r ~ stim,
            gabor)
plot(m6)

## If you do not want to see how the model fits at the individual level you can specify a subset of
## predictors as variables of interest
plot(m6, vs = c('duration'))

## You can get some information about the model parametrization by looking at the model object
m6

####################################################################################################
## meta-d'

## You can try to fit the meta-d' generalization like so:
m7 = bhsdtr(c(metad ~ duration * order + (duration | id), thr ~ order + (1 | id)),
            r ~ stim,
            gabor)
plot(m7)
## This model also seems to fit quite well. It has to fit at least as well as the regular SDT
## version regardless if it is true or if it makes sense.

## The metad parameter is two-dimensional. The first dimension is d', the second dimension is
## meta-d'.
samples(m7, 'metad')

## Let's say you want to obtain the jmap point estimate of the meta-d'/d' ratio for the 32 ms x
## DECISION-RATING condition. All you have to do is this:
s = samples(m7, 'metad')
s[,2,'32 ms:DECISION-RATING'] / s[,1,'32 ms:DECISION-RATING']
## If there were posterior samples in the samples object you would obtain samples from the posterior
## meta-d'/d' distribution.

## If I was forced to do this, this is how I would approach the problem of making inference about
## the meta-d' effects or the meta-d'/d' effects in this case:

## To make the chains go faster we will simplify the problem
m8 = bhsdtr(c(metad ~ order + (order | id), thr ~ order + (order | id)),
            r ~ stim,
            gabor[gabor$duration == '32 ms',], method = 'stan',
            chains = 7)
## There are problems with this model (>100 divergent transitions) but we don't have time to
## fine-tune it (the main thing to try is altering the priors and increasing the adapt_delta
## parameter).
## 
## saveRDS(m8, file = 'E:/assc25m8.rds')
m8 = readRDS('E:/assc25m8.rds')

s = samples(m8, 'metad')
## Let's look for a reason to introduce the meta-d' parameter
apply(s[,2,] / s[,1,], 2, function(x)quantile(x, c(.025, .975)))
##       DECISION-RATING RATING-DECISION
## 2.5%        0.2331649       0.3283728
## 97.5%       0.7862132       1.5515857
##
## The meta-d' / d' ratio seems to be different from 1 in the DECISION-RATING condition

####################################################################################################
## Average Threshold Spread (ATS) as a measure of metacognitive judgement bias

m9 = bhsdtr(c(dprim ~ 1 + (1 | id), thr ~ 1 + (1 | id)),
            r ~ stim,
            gabor[gabor$duration == '32 ms' & gabor$order == 'DECISION-RATING',],
            method = 'stan',
            chains = 7)
## saveRDS(m9, file = 'E:/assc25m9.rds')
m9 = readRDS('E:/assc25m9.rds')

## This is how you can obtain posterior samples of threshold for each person (or level of some other grouping factor)
s = samples(m9, 'thr', group = 1)
## We want separate posterior samples of ATS for each condition x person
ats = matrix(ncol = dim(s)[3], nrow = dim(s)[1])
## and we would like to preserve the names of the conditions x persons - this makes it easier to
## match the posterior samples to individual participants
colnames(ats) = dimnames(s)[[3]]
for(i in 1:(dim(s)[3]))
    ## The fourth threshold is the middle one in this case, since there are 7 of them
    ats[,i] = apply(s[,,i], 1, function(x)mean(abs(x - x[4])))

mean.rating = aggregate(rating ~ id, gabor[gabor$duration == '32 ms' & gabor$order == 'DECISION-RATING',], mean)
mean.ats = apply(ats, 2, mean)
all(names(mean.ats) == mean.rating$id)
## The names match
plot(mean.rating$rating ~ mean.ats)
cor(mean.rating$rating, mean.ats)
## -0.8123227

## ATS corresponds closely to average ratings but, as long as the model is true, it does not depend
## on type 1 sensitivity
