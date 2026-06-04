library(shiny)
library(ggplot2)
library(SciViews)
library(scales)
E<-100
pr<-0.1


# Define the functions
maxAUC <- function(FOM, pr) {
  p<-E*pr
  q<-E-p
  a <- FOM * p
  if (FOM >= pr) {
    result <- (2 * p * q - (p^2 / a - p) * (p - a)) / (2 * p * q)
  } else {
    result <- (2 * p * q - (E - p) * (p - a)) / (2 * p * q)
  }
  return(result)
}

minAUC <- function(FOM, pr) {
  result <- (1 + FOM) / 2
  result <- max(pr,result)
  return(result)
}

maxAFOM<-function(FOM, pr){
  p<-E*pr
  q<-E-p
  Tot<-p*(1/2-ln(pr)-pr/2)
  
  numerator<- -p*(pr/2-ln(FOM/pr)-1+FOM^2/2)
  # result<-max(numerator/Tot,pr)
  result<-numerator/Tot
  return(result)
}
minAFOM<-function(FOM, pr){
  p<-E*pr
  q<-E-p
  h<-function(x){pmin(x, (p-FOM*p)/(E-FOM*p)*(x-FOM*p)+FOM*p)}
  f<-function(x){h(x)/(x-h(x)+p)}
  numerator<-integrate(f, lower=0, upper=E)$value-p*pr/2
  Tot<-p*(1/2-ln(pr)-pr/2)
  return(numerator/Tot)
}

# Define UI for the app
ui <- fluidPage(
  titlePanel("Interactive Plot with AUC Functions"),
  
  sidebarLayout(
    sidebarPanel(
      sliderInput("pr", "Prevalence (pr):", min = 0, max = 1, value = 0.1),
    ),
    
    mainPanel(
      plotOutput("aucPlot"),
      plotOutput("afomPlot")
    )
  )
)

# Define server logic
server <- function(input, output) {
  
  # Reactive expression to calculate the AUC values based on inputs
  output$aucPlot <- renderPlot({
    FOM_values <- seq(input$pr, 1, length.out = 100)
    
    # Get values for maxAUC and minAUC
    max_auc_values <- sapply(FOM_values, function(FOM) maxAUC(FOM, input$pr))
    min_auc_values <- sapply(FOM_values, function(FOM) minAUC(FOM, input$pr))
    
    # Create a data frame for plotting
    auc_data <- data.frame(
      FOM = rep(FOM_values, 2),
      AUC = c(max_auc_values, min_auc_values),
      Type = rep(c("maxAUC", "minAUC"), each = length(FOM_values))
    )
    
    # Plot using ggplot2
    ggplot(auc_data, aes(x = FOM, y = AUC, color = Type)) +
      geom_line(size = 1) +
      geom_hline(yintercept = 0.9, linetype="dashed")+
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
      labs(title = "AUC vs FOM", x = "FOM", y = "AUC") +
      scale_x_continuous(limits = c(0, 1), expand = c(0, 0), breaks = seq(0, 1, by = 0.2), labels = percent_format()) +
      scale_y_continuous(limits = c(0, 1), expand = c(0, 0), breaks = seq(0, 1, by = 0.2), labels = percent_format()) + 
      theme_minimal() +
      theme(legend.title = element_blank())+
      theme(aspect.ratio=1)
  })
  output$afomPlot <- renderPlot({
    FOM_values <- seq(input$pr, 1, length.out = 100)
    
    # Get values for maxAUC and minAUC
    max_auc_values <- sapply(FOM_values, function(FOM) maxAFOM(FOM, input$pr))
    min_auc_values <- sapply(FOM_values, function(FOM) minAFOM(FOM, input$pr))
    
    # Create a data frame for plotting
    auc_data <- data.frame(
      FOM = rep(FOM_values, 2),
      AUC = c(max_auc_values, min_auc_values),
      Type = rep(c("maxAUC", "minAUC"), each = length(FOM_values))
    )
    
    # Plot using ggplot2
    ggplot(auc_data, aes(x = FOM, y = AUC, color = Type)) +
      geom_line(size = 1) +
      geom_hline(yintercept = 0.9, linetype="dashed")+
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
      labs(title = "AUC vs FOM", x = "FOM", y = "AUC") +
      scale_x_continuous(limits = c(0, 1), expand = c(0, 0), breaks = seq(0, 1, by = 0.2), labels = percent_format()) +
      scale_y_continuous(limits = c(0, 1), expand = c(0, 0), breaks = seq(0, 1, by = 0.2), labels = percent_format()) + 
      theme_minimal() +
      theme(legend.title = element_blank())+
      theme(aspect.ratio=1)
  })
}

# Run the application
shinyApp(ui = ui, server = server)

