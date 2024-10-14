using Distributions
using Expectations

xs
ps


d = MultivariateDiscreteNonParametric(xs, ps)
E = expectation(d)