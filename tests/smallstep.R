library(calibrar)
library(optimx)
path = NULL # NULL to use the current directory
ARPM = calibrar_demo(path=path, model="PoissonMixedModel", L=5, T=100) 
setup = calibration_setup(file=ARPM$setup)
observed = calibration_data(setup=setup, path=ARPM$path)
forcing = as.matrix(read.csv(file.path(ARPM$path, "master", "environment.csv"), row.names=1))

# real parameters
real = ARPM$par

# fixing some parameters
ARPM$guess$sd = NULL
ARPM$guess$mu_ini = NULL
ARPM$lower$sd = NULL
ARPM$lower$mu_ini = NULL
ARPM$upper$sd = NULL
ARPM$upper$mu_ini = NULL

skeleton = as.relistable(ARPM$guess)

run_model = function(par, forcing) {
  par = relist(flesh=par, skeleton = skeleton)
  par$sd = real$sd
  par$mu_ini = real$mu_ini
  output = calibrar:::.PoissonMixedModel(par=par, forcing=forcing)
  output = c(output, list(gammas=par$gamma)) # adding gamma parameters for penalties
  return(output)
}

# objective function
fn = calibration_objFn(model=run_model, setup=setup, observed=observed, forcing=forcing, aggregate=TRUE)
gr = function(par) gradient(fn=fn, x=par, method="forward")

guess = unlist(ARPM$guess)
lower = unlist(ARPM$lower)
upper = unlist(ARPM$upper)

# message printing
msg = function(test, deltat=NULL) {
  if(is.null(deltat)) deltat = test$elapsed
  cat(sprintf("value=%0.3f | time=%s | counts=(%s)\n", 
              test$value, format(deltat, units="seconds", digits=3), 
              paste(test$counts, collapse=",")))
  return(invisible())
}

maxit = 2000
maxfeval = 20000 # There's a typo in the documentation, it's not 'maxfevals'.

cat(sprintf("The minimum value of the function is %0.3f\n", ARPM$value))

# smallstep tests ---------------------------------------------------------

# original Rvmmin
test0 = Rvmmin.test(par=guess, fn=fn, gr=gr, lower=lower, upper=upper, 
               control=list(maxit=maxit, maxfeval=maxfeval))
msg(test0)

# smallstep = .Machine$double.xmin
test1 = Rvmmin.test(par=guess, fn=fn, gr=gr, lower=lower, upper=upper, 
                            control=list(maxit=maxit, maxfeval=maxfeval, 
                                         smallstep=.Machine$double.xmin))
msg(test1)

# not saving unfeasible points
test2 = Rvmmin.test(par=guess, fn=fn, gr=gr, lower=lower, upper=upper, 
                            control=list(maxit=maxit, maxfeval=maxfeval, nosave=TRUE))
msg(test2)

# smallstep = .Machine$double.xmin, not saving unfeasible points
test3 = Rvmmin.test(par=guess, fn=fn, gr=gr, lower=lower, upper=upper, 
                            control=list(maxit=maxit, maxfeval=maxfeval, smallstep=.Machine$double.xmin,
                                         nosave=TRUE))
msg(test3)


# Internal numerical gradient tests ---------------------------------------

# here we use internal numerical gradient computation
# original Rvmmin (with internal gradient)
test0b = Rvmmin.test(par=guess, fn=fn, lower=lower, upper=upper, 
                             control=list(maxit=maxit, maxfeval=maxfeval))
msg(test0b)
# smallstep = .Machine$double.xmin (with internal gradient)
test1b = Rvmmin.test(par=guess, fn=fn, lower=lower, upper=upper, 
                             control=list(maxit=maxit, maxfeval=maxfeval, 
                                          smallstep=.Machine$double.xmin))
msg(test1b)



