### Implementation of a toy GUTS model in order to play
##  with the parameters and have a more practical
##  understanding of the functioning

library(deSolve)
library(ggplot2)
library(reshape)
library(cowplot)
library(plotly)

# function to be passed to the solver
# it contains all the expressions
guts_equations_sd <- function(t, state, parms, approxfun){
  with(as.list(c(state, parms)),{
    inputCw <- approxfun(x = conc.time, y = conc.Cw, rule = 2)
    cw <- cbind(inputCw(t)) # external concentration; mg/L
    dCi <- ke * (Kiw * cw - Ci)
    dDi <- kr * (Ci - Di)
    hz <- bi * max(0, Di - mi) + hb
    dS <- -hz * S
    
    list(c(dCi, dDi, dS))
  })
}


# implementation for the IT case more difficult, need to find
# right way to import the loglogistic function in R
# TO BE DONE
guts_equations_it <- function(t, state, parms, approxfun){
  with(as.list(c(state, parms)),{
    inputCw <- approxfun(x = conc.time, y = conc.Cw, rule = 2)
    cw <- cbind(inputCw(t))
    dCi <- ke * (Kiw * cw - Ci)
    dDi <- kr * (Ci - Di)
    dS <- 0
    list(c(dCi, dDi, dS))
  })
}


compute_survival <-function(time, state, parms, ode_res, approxfun){
  damagetemp<-0
  output <- ode_res
  with(as.list(c(state, parms)),{
    inputD <- approxfun(x=ode_res$time, ode_res$Di, rule=2)
    #damage=cbind(inputD(time))
    for (i in 1:length(time)){
    Sb <- exp(-hb*time[i])
    damage<-inputD(time[i])
    if (damage>damagetemp){damagetemp<-damage}
    Surv <- Sb * (1 - cumulative_loglogistic(damagetemp, mi, beta))
    output$S[i] <- Surv
    }
    return(as.data.frame(output))
  })
}

cumulative_loglogistic <- function(x, m, beta){
  cumulative <- 1. / (1 + (x / m)^(-beta))
  return(cumulative)
}


# Parameters:
# ke      = elimination rate constant. how rapidly the concentrations inside
#           and outside even
# Kiw     = bioconcentration factor, it measures affinity of body tissues wrt 
#           affinity with the external medium 
#           (in steady state Body residual = Kiw * external conc.)
# kr      = damage repair rate. Rate at which damage reaches equilibrium with 
#           body residual
# mi/mw   = median of threshold distribution (for IT case)
# beta/Fs = spread factor of the threshold distribution
# bi/bw   = killing rate. Determines the steepness of the probability to die
#           when going above threshold
# hb      = background hazard rate


# State variables:
# Ci = internal concentration
# Di = damage value
# S  = survival probability


# System parameters
# Cw = external concentration
# N  = number of organisms ?


parms <- c(ke=0.9, 
           Kiw=1., 
           kr=0.3,
           k=0.3,
           mi=2.5,
           beta=2,
           bi=0.1,
           hb=0)

state_vars <- c(Ci=0,
                Di=0,
                S=1)
state_vars2 <- c(Di=0,
                S=1)


time <- seq(0,7,by=0.1)
extc <- 2  # external concentration
Cw <- c(rep(0,length(time)))
Cw[c(1:40)]<-5
Cw[c(40:70)]<-0
conce <- data.frame(time=time, Cw=Cw)

guts_red_it <- function(t, state, parms, approxfun){
  ## as it is a toy model, the arguments are called this way, keeping names
  with(as.list(c(state, parms)),{
    inputCw <- approxfun(x = conc.time, y = conc.Cw, rule = 2)
    cw <- cbind(inputCw(t))
    dDi <- k * (cw - Di)
    dS <- 0
    list(c(dDi, dS))
  })
}

solve <- function(parameters, initial, time,modeltype="SD"){
  if (modeltype=="SD"){
    out_sd <- ode(y=initial, times=time, func=guts_equations_sd, 
                  parms = parameters,
                  approxfun = approxfun)
    out <- as.data.frame(out_sd)
  }
  if (modeltype == "IT"){
    out_it <- ode(y=state_vars, times=time, func=guts_equations_it, 
                  parms = c(parms, conc=conce),
                  approxfun = approxfun)
    out_it <- as.data.frame(out_it)
    out <- compute_survival(time, initial, parms, out_it, approxfun)
  }
  return(out)
}

solve_and_plot <- function(parameters, initial, time, modeltype="SD"){
  plotpdf<-NULL
  if (modeltype=="SD"){
      out_sd <- ode(y=initial, times=time, func=guts_equations_sd, 
                parms = parameters,
                approxfun = approxfun)
      out <- as.data.frame(out_sd)
  }
  if (modeltype == "IT"){
    out_it <- ode(y=state_vars, times=time, func=guts_equations_it, 
                  parms = c(parms, conc=conce),
                  approxfun = approxfun)
    out_it <- as.data.frame(out_it)
    out <- compute_survival(time, state_vars, parms, out_it, approxfun)
    xvals <- seq(0,log10(parms["mi"]*10^(3*parms["mi"])),by=0.1)
    yvals <- actuar::dllogis(xvals,shape = parameters$beta, scale = parameters$mi)
    pdfdata<-as.data.frame(cbind(xvals, yvals))
    plotpdf <- ggplot(pdfdata) + geom_line(aes(x=xvals, y=yvals), size=1) +
      xlab("Threshold") + ylab("pdf")
  }
  plotCi <- ggplot(out) + geom_line(aes(x=time, y=Ci), size=1) +
    xlab("Time [a.u.]") + ylab("Internal concentration [a.u.]")
  plotDi <- ggplot(out) + geom_line(aes(x=time, y=Di), size=1) +
    xlab("Time [a.u.]") + ylab("Internal Damage [a.u.]")
  plotS <- ggplot(out) + geom_line(aes(x=time, y=S), size=1) +
    xlab("Time [a.u.]") + ylab("Survival probability")
  if (is.null(plotpdf)){
    p<-plot_grid(plotCi, plotDi, plotS)
  }
  else {
    p<-plot_grid(plotCi, plotDi, plotS, plotpdf)
    }
  return(p)
}


solve_and_plot2 <- function(parameters, initial, time, modeltype="IT"){
  plotpdf<-NULL
  if (modeltype=="SD"){
    out_sd <- ode(y=initial, times=time, func=guts_equations_sd, 
                  parms = parameters,
                  approxfun = approxfun)
    out <- as.data.frame(out_sd)
  }
  if (modeltype == "IT"){
    out_it <- ode(y=state_vars2, times=time, func=guts_red_it, 
                  parms = c(parms, conc=conce),
                  approxfun = approxfun)
    out_it <- as.data.frame(out_it)
    out <- compute_survival(time, state_vars2, parms, out_it, approxfun)
    xvals <- seq(0,log10(parms["mi"]*10^(3*parms["mi"])),by=0.1)
    yvals <- actuar::dllogis(xvals,shape = parameters$beta, scale = parameters$mi)
    pdfdata<-as.data.frame(cbind(xvals, yvals))
    plotpdf <- ggplot(pdfdata) + geom_line(aes(x=xvals, y=yvals), size=1) +
      xlab("Threshold") + ylab("pdf")
  }
  plotDi <- ggplot(out) + geom_line(aes(x=time, y=Di), size=1) +
    xlab("Time [a.u.]") + ylab("Internal Damage [a.u.]")
  plotS <- ggplot(out) + geom_line(aes(x=time, y=S), size=1) +
    xlab("Time [a.u.]") + ylab("Survival probability")
  if (is.null(plotpdf)){
    p<-plot_grid(plotCi, plotDi, plotS)
  }
  else {
    p<-plot_grid(plotDi, plotS, plotpdf)
  }
  return(p)
}

solve_and_plot(c(parms, conc=conce), state_vars, time, modeltype = "SD")
solve_and_plot(c(parms, conc=conce), state_vars, time, modeltype = "IT")
solve_and_plot2(c(parms, conc=conce), state_vars, time, modeltype = "IT")

