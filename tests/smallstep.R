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

control = list(maxit=20000, eps=sqrt(.Machine$double.eps))

# objective function
fn = calibration_objFn(model=run_model, setup=setup, observed=observed, forcing=forcing, aggregate=TRUE)
gr = function(par) gradient(fn=fn, x=par, method="forward")

guess = unlist(ARPM$guess)
lower = unlist(ARPM$lower)
upper = unlist(ARPM$upper)

# original Rvmmin
test0 = optimx::Rvmmin(par=guess, fn=fn, gr=gr, lower=lower, upper=upper, 
               control=list(maxit=20000, eps=sqrt(.Machine$double.eps)))


ori = Rvmmin_ori(par=guess, fn=fn, lower=lower, upper=upper, control=control)
ori2 = Rvmmin_ori(par=guess, fn=fn, gr=gr, lower=lower, upper=upper, control=control)
mod = Rvmmin(par=guess, fn=fn, lower=lower, upper=upper, control=control)
mod2 = Rvmmin(par=guess, fn=fn, gr=gr, lower=lower, upper=upper, control=control)
test = optim2(par=guess, fn=fn, lower=lower, upper=upper, control=control, method="Rvmmin")


