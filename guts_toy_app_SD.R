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

# function to be passed to deSolve::ode
# it contains all the expressions
# Solves the model with a fixed threshold
# parameters:
#  - t : time
#  - state: initial conditions of the state variables
#  - parms: parameters (including forcing ones)
#  - approxfun: auxiliary function to allow time dependent ext. concentration

guts_red_sd <- function(t, state, parms, approxfun){
  ## as it is a toy model, the arguments are called this way, keeping names
  with(as.list(c(state, parms)),{
    inputCw <- approxfun(x = conc.time, y = conc.Cw, rule = 2)
    cw <- cbind(inputCw(t))
    dDi <- k * (cw - Di)
    hz <- bi * max(0, Di - mi) + hb
    dS <- -hz * S
    
    list(c(dDi, dS))
  })
}

# Solve the differential equation for GUTS SD or IT
solve <- function(parameters, initial, time){
    out_sd <- ode(y=initial, times=time, func=guts_red_sd,
                  parms = parameters,
                  approxfun = approxfun)
    out <- as.data.frame(out_sd)
  return(out)
}

## MAIN CODE for the application
# Base code taken from internet snippets
# call the user interface object
ui <- dashboardPage(
  dashboardHeader(title="GUTS SD model"),
  dashboardSidebar(sliderInput("Cw","External Concentration", min=0, max=10, step=0.1, value=5),
                   sliderInput("k","Damage Repair/Elimination rate", min=0., max=5, step=0.1, value=0.5),
                   sliderInput("mi","Threshold", min=0., max=10, step=0.1, value=3),
                   sliderInput("bi","Killing rate", min=0., max=2, step=0.1, value=0.5),
                   sliderInput("hb","Background hazard", min=0., max=0.01, step=0.001, value=0.001)
                   ),
  dashboardBody(
    fluidRow(box(plotOutput("guts"), height = 800, width = 750))
  ))

# Function that actually call the computations and produces the plots
server <- function(input, output, session) {
  output$guts <- renderPlot({
    parms <- c(input$k, input$mi, input$bi, input$hb)
    names(parms)<-c("k","mi","bi","hb")
    time <- seq(0,14,by=0.01)
    Cw <- c(rep(0,length(time)))
    Cw[c(1:700)]<-input$Cw
    conc <- data.frame(time=time, Cw=Cw)
    state_vars <- c(Di=0,S=1)
    out_sd <- solve(c(parms,conc=conc), state_vars, time)
    plotCi <- ggplot(out_sd) + geom_line(aes(x=time, y=Cw), linewidth=1) +
      xlab("Time [a.u.]") + ylab("External concentration [a.u.]")+
      theme(text=element_text(size=14))
    plotDi <- ggplot(out_sd,aes(x=time, y=Di)) + geom_line(linewidth=1) +
      xlab("Time [a.u.]") + ylab("Internal Damage [a.u.]")+
      geom_hline(yintercept=input$mi, linetype=2)+
      theme(text=element_text(size=14))
    plotS <- ggplot(out_sd) + geom_line(aes(x=time, y=S), linewidth=1) +
      xlab("Time [a.u.]") + ylab("Survival probability")+ylim(0,1)+
      theme(text=element_text(size=14))
    plot_grid(plotCi, plotDi, plotS, ncol = 1)
  }, height = 750)
}

# launch the applet
shinyApp(ui, server)
