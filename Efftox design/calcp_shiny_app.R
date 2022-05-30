library(shiny)
library(tidyverse)
library(graphics)
# rmarkdown::render_delayed(
shinyApp(
  ui <- shinyUI(fluidPage(
    withMathJax(), 
    titlePanel("Eff-Tox Trade-Off Contour"), 
    sidebarLayout(
      sidebarPanel (
        numericInput(inputId = 'pistar1E',label = '\\((\\pi^*_{1,E},0)\\): Input the value of
                     \\(\\pi^*_{1,E}\\)',value=0.5,min=0,max=1, step=0.00001),
        numericInput(inputId = 'pistar2T',label = '\\((1,\\pi^*_{2,T})\\): Input the value of
                     \\(\\pi^*_{2,T}\\)',value=0.65,min=0,max=1,step=0.00001),
        numericInput(inputId = 'pistar3E',label = '\\((\\pi^*_{3,E},\\pi^*_{3,T})\\): Input the value of
                   \\(\\pi^*_{3,E}\\)',value=0.7,min=0,max=1,step=0.00001),
        numericInput(inputId = 'pistar3T',label = '\\((\\pi^*_{3,E},\\pi^*_{3,T})\\): Input the value of
                     \\(\\pi^*_{3,T}\\)',value=0.25,min=0,max=1,step=0.00001),
        sliderInput(inputId = "utility_lower",
                    label = "Utility_lower",
                    min=-3,
                    max=0,
                    value = -3,
                    step = 0.1),
        sliderInput(inputId = "utility_upper",
                    label = "Utility_upper",
                    min=0,
                    max=3,
                    value = 3,
                    step = 0.1)        
        
      ),
      mainPanel(
        textOutput("pvalue")
        ,plotOutput("contour",width ="100%")
        ,plotOutput("contour2",width ="100%")
        ,textOutput("test")
      )
    )
  )
  ),
  server <- function(input, output) {
    pi1E=reactive({input$pistar1E})
    pi2T=reactive({input$pistar2T})
    pi3E=reactive({input$pistar3E})
    pi3T=reactive({input$pistar3T})
    utility_lower=reactive({input$utility_lower})
    utility_upper=reactive({input$utility_upper})
    pcalc= function(p, pi1E, pi2T, pi3E, pi3T) {
      a <- ((1 - pi3E) / (1 - pi1E))
      b <- pi3T / pi2T
      return(a^p + b^p - 1)
    }
    
    p=reactive({stats::uniroot(pcalc, interval = c(0, 100),pi1E=pi1E(),pi2T=pi2T(),pi3E=pi3E(),pi3T=pi3T())})
    pv<-reactive({p()$root})
    
    
    efftox_get_tox <- function(eff, util, p, eff0, tox1) {
      a = ((1 - eff) / (1 - eff0))
      return(tox1 * ((1 - util)^p - a^p)^(1 / p))
    }
    
    
    # utival=seq(-1,1,by=0.1)
    utival=reactive({seq(utility_lower(),utility_upper(),0.1)})
    effgrid<-seq(0,1,length.out=1000)
    
    toxval=reactive({sapply(utival(),function(u) efftox_get_tox(eff = effgrid,util = u,p=pv(),eff0=pi1E(),tox1=pi2T()))})
    
    dataf=reactive({data.frame(effgrid=rep(effgrid,times=length(utival())),
                               toxval=as.numeric(toxval()),
                               utival=rep(utival(),each=length(effgrid)))
    })
    
    tox_vals_fix = reactive({efftox_get_tox(eff = effgrid,util = 0, pv(), eff0=pi1E(),tox1=pi2T())})
    dataf2 = reactive({data.frame(effgrid=rep(effgrid,times=length(utival())),
                                  toxval=rep(as.numeric(tox_vals_fix(),each=length(utival()))),
                                  utivals = rep(0,1000,each=length(utival())))})
    
    # plt=reactive({ggplot(data = dataf(),aes(x=effgrid,y=toxval,group=as.factor(utival)))+geom_line(size = 0.5, alpha = 0.25,col="red")+
    #     xlim(0,1)+ylim(0,1)+xlab('Prob(Efficacy)')+ylab('Prob(Toxicity)')+geom_line(data = dataf2(), size = 1,col="red")})
    
    output$contour<-renderPlot({
      ggplot(data = dataf(),aes(x=effgrid,y=toxval,group=as.factor(utival)))+geom_line(size = 1, alpha = 0.25)+
        xlim(0,1)+ylim(0,1)+xlab('Prob(Efficacy)')+ylab('Prob(Toxicity)')+
        geom_line(aes(y=efftox_get_tox(effgrid,0,pv(),eff0=pi1E(),tox1=pi2T()),x=effgrid),lwd=1.5,col="red")
      # +
      #   geom_point(aes(x=pi1E(),y=0),col="blue")+
      #   geom_point(aes(x=1,y=pi2T()),col="blue")+
      #   geom_point(aes(x=pi3E(),y=pi3T()),col="blue")
    })
    
    
    # output$test<-renderPrint({
    #   length(utival1())
    # })
    
    x1=seq(-1,0,0.01)
    xlab1=seq(0,1,0.01)
    y1=seq(0,1,0.01)
    
    
    z1=reactive({1-(outer(abs(x1)^pv(),y1^pv(),"+"))^(1/pv())})
    
    output$pvalue<-renderText(paste0("The value of p: ", round(pv(),5)))
    
    output$contour2<-renderPlot({
      # plot(NULL,ylim=c(0,1),xlim=c(0,1),ylab="P(T)",xlab="P(E)")
      plt<-contour(xlab1,y1,z1(),labcex=1.5,lwd=2)
      plt
      # for (u in utival){
      # ylab=efftox_get_tox(xlab,u,pv(),eff0=pi1E(),tox1=pi2T())
      # points(y=ylab,x=xlab,lwd=2,type="l")
      # }
      lines(y=efftox_get_tox(xlab1,0,pv(),eff0=pi1E(),tox1=pi2T()),x=xlab1,lwd=2,col="red")
    })
    
    # pass the value of p to global environment for further use in stan
    observe({
      assign(
        x="pval",
        value = pv(),
        envir = .GlobalEnv
      )
    })
    observe({
      assign(
        x="input",
        value = c(pi1E=pi1E(),pi2T=pi2T(),pi3E=pi3E(),pi3T=pi3T()),
        envir = .GlobalEnv
      )
    })
    observe({
      assign(
        x="utilitygrid",
        value = utival(),
        envir = .GlobalEnv
      )
    })
  }
)