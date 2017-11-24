# Issues and Future Directions in `cnqr` Package

* Integration of joint density is slow when the most recent vine was at least 2 iterations ago.
	* Example: Integrate-out 2 variables to get joint density `fX` -- evaluates in 25 seconds. Get conditional distribution by computing two integrals of `fX` (one up to the observation, one over the whole range) -- about 150 evals per integral, means 150x2x25s = ~2h per observation. With 100 observations, that's 8.7 days.
	* Fix 1: When integrating a variable in `fX`, the denominator integral goes from 0 to 1 -- which is the marginal distribution _without_ the variable we're integrating out. That's used in the next sequence if the next subset still doesn't form a vine.
	* Fix 2: The integrals themselves might be sped up by some sort of Gauss-Legendre integration.
	* Entirely different idea: Why not determine the order of predictors first, then fit a vine for each set of predictors? This way, the most we'd have to integrate is once each time. After all, we're only looking for sequential estimates of one predictor given the previous, and however it's done should be fine -- we're _not_ necessarily looking to model the joint distribution of the predictors. But, we'd be fitting that many more parameters, which is likely a problem if the goal is estimation (what are the problems?), but do those problems still hold if we're ultimately forecasting?

* Bayesian Network order: currently, `lm()` is used to get conditional values for which correlation is computed. This is crude, because:
	* the mean line is probably not linear. If I really want to optimize the order of sequential partial correlations, then I should probably fit copulas and get the kendall's tau for the conditional probability transforms.
	* I don't care about overall correlation. The analogue of this that I would care about is the predictive power a predictor has on the upper part of `Y` (of course, both conditioned on previous predictors).
	* But, possibly what I really care about is achieving an order where the simplifying assumption is the least far from the truth -- because then it follows that a copula really does link the conditional response with the conditional predictors.

* Should include estimates of covariance matrix of cnqr estimator.

* When specifying weights across quantile levels, `scoreq()`'s `wtau` argument should really also accept a vector corresponding to the inputted taus, instead of only accepting a function.

* Add bibtex citation to the package

* Tiffany was telling me about Ropensci for reviewing the package, and getting a "stamp of approval". Two people will review the package to give feedback. Need `testthat` tests for this. It's a nice stamp of approval. In addition, Journal of Open Source Software is another option besides Journal of Statistical Software, which I can submit for review by checking off a box during Ropensci application (just requires a small writeup). Good to add to publication list. 

* Replace the name `scorer` with the name `psr` for "proper scoring rule". Perhaps to actually evaluate scores, one can type `psr_eval` instead of the original `scoreq` or even the later `score_eval`.

* `cnqr` function has `cdf` argument that doesn't make much sense. I think the documentation means to say that the cdfs are only needed for the predictors? I think the purpose of this is to make the `cnqr` function do more things, particularly do the PIT transformation within it. But do I need the cdf of Y? I don't think so, yet at the same time, I seem to vaguely remember needing to include it for significant improvements in speed. 