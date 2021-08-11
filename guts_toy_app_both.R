### Implementation of a toy GUTS model in order to play
##  with the parameters and have a more practical
##  understanding of the functioning
##
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
#
# State variables:
# Ci = internal concentration
# Di = damage value
# S  = survival probability
#
# System parameters
# Cw = external concentration
# N  = number of organisms ?

### TODO:
# Optimize the design of the sliders with more intelligent parameters



library(deSolve)
library(ggplot2)
library(reshape)
library(cowplot)
library(plotly)
library(shiny)
library(shinydashboard)
library(shinyWidgets)

# function to be passed to deSolve::ode
# it contains all the expressions
# Solves the model with a fixed threshold
# parameters:
#  - t : time
#  - state: initial conditions of the state variables
#  - parms: parameters (including forcing ones)
#  - approxfun: auxiliary function to allow time dependent ext. concentration
guts_equations_sd <- function(t, state, parms, approxfun){
  ## as it is a toy model, the arguments are called this way, keeping names
  with(as.list(c(state, parms)),{
    inputCw <- approxfun(x = conc.time, y = conc.Cw, rule = 2)
    cw <- cbind(inputCw(t)) 
    dCi <- ke * (Kiw * cw - Ci)
    dDi <- kr * (Ci - Di)
    hz <- bi * max(0, Di - mi) + hb
    dS <- -hz * S
    
    list(c(dCi, dDi, dS))
  })
}


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
}


# Solve the differential equation for GUTS SD or IT
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
    out <- compute_survival(time, state_vars, parms, out_it, approxfun)
  }
  return(out)
}

## MAIN CODE for the application
# Base code taken from internet snippets
# call the user interface object
ui <- dashboardPage(
  dashboardHeader(title="GUTS SD and IT models"),
  dashboardSidebar(box(sliderTextInput(inputId = "model", label = "Model Type",
                                       choices=c("SD", "IT"), selected="SD"),
                       title = "Model type",
                       width = 150, height = 150),
                   conditionalPanel(condition = "input.model==SD", 
                                    selectInput(sliderInput("Cw","External Concentration", min=0, max=10, step=0.1, value=2),
                                                sliderInput("ke","Elimination rate", min=0., max=5, step=0.1, value=0.5),
                                                sliderInput("Kiw","Bioconcentration factor", min=0., max=2, step=0.1, value=1),
                                                sliderInput("kr","Damage repair rate", min=0., max=5, step=0.1, value=0.5),
                                                sliderInput("mi","Threshold", min=0., max=10, step=0.1, value=1.5),
                                                sliderInput("bi","Killing rate", min=0., max=2, step=0.1, value=0.5),
                                                sliderInput("hb","Background hazard", min=0., max=0.01, step=0.001, value=0.001))),
                   conditionalPanel(condition = "input.model==IT", 
                                    selectInput(sliderInput("Cw","External Concentration", min=0, max=10, step=0.1, value=2),
                                                sliderInput("ke","Elimination rate", min=0., max=5, step=0.1, value=0.5),
                                                sliderInput("Kiw","Bioconcentration factor", min=0., max=2, step=0.1, value=1),
                                                sliderInput("kr","Damage repair rate", min=0., max=5, step=0.1, value=0.5),
                                                sliderInput("mi","Threshold", min=0., max=10, step=0.1, value=1.5),
                                                sliderInput("beta","Shape parameter", min=0., max=5, step=0.01, value=0.5),
                                                sliderInput("hb","Background hazard", min=0., max=0.01, step=0.001, value=0.001)))
                   ),
  dashboardBody(
    fluidRow(box(plotOutput("guts"), height = 800, width = 750))
  ))

# Function that actually call the computations and produces the plots
server <- function(input, output, session) { 
  # need to use "reactive" to process the condition on the ui
  output$guts <- renderPlot({
    parms <- c(input$ke, input$Kiw, input$kr, input$mi, input$bi, input$hb)
    names(parms)<-c("ke","Kiw","kr","mi","bi","hb")
    time <- seq(0,100,by=0.1)
    Cw <- c(rep(0,length(time)))
    Cw[c(1:200)]<-input$Cw
    conc <- data.frame(time=time, Cw=Cw)
    state_vars <- c(Ci=0,Di=0,S=1)
    out_sd <- solve(c(parms,conc=conc), state_vars, time)
    plotCi <- ggplot(out_sd) + geom_line(aes(x=time, y=Ci), size=1) +
      xlab("Time [a.u.]") + ylab("Internal concentration [a.u.]")
    plotDi <- ggplot(out_sd) + geom_line(aes(x=time, y=Di), size=1) +
      xlab("Time [a.u.]") + ylab("Internal Damage [a.u.]")
    plotS <- ggplot(out_sd) + geom_line(aes(x=time, y=S), size=1) +
      xlab("Time [a.u.]") + ylab("Survival probability")
    plot_grid(plotCi, plotDi, plotS, ncol = 1)
  }, height = 750)
}

# launch the applet
shinyApp(ui, server)
