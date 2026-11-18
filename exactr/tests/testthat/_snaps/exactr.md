# binom.test works

    Code
      binom.test(c(800, 10))
    Output
      
      	Exact binomial test
      
      data:  c(800, 10)
      number of successes = 800, number of trials = 810, p-value < 2.2e-16
      alternative hypothesis: true probability of success is not equal to 0.5
      95 percent confidence interval:
       0.9774134 0.9940643
      sample estimates:
      probability of success 
                   0.9876543 
      

# poisson.test works

    Code
      poisson.test(137, 24.19893)
    Output
      
      	Exact Poisson test
      
      data:  137 time base: 24.19893
      number of events = 137, time base = 24.199, p-value < 2.2e-16
      alternative hypothesis: true event rate is not equal to 1
      95 percent confidence interval:
       4.753125 6.692709
      sample estimates:
      event rate 
        5.661407 
      
    Code
      poisson.test(c(11, 6 + 8 + 7), c(800, 1083 + 1050 + 878))
    Output
      
      	Comparison of Poisson rates
      
      data:  c(11, 6 + 8 + 7) time base: c(800, 1083 + 1050 + 878)
      count1 = 11, expected count1 = 6.7174, p-value = 0.07967
      alternative hypothesis: true rate ratio is not equal to 1
      95 percent confidence interval:
       0.8584264 4.2772659
      sample estimates:
      rate ratio 
        1.971488 
      

