library(shiny)
library(ggplot2)
library(gridExtra)

# Define la interfaz de usuario
ui <- fluidPage(
  titlePanel("Comparación de Matrices con Raster"),
  
  sidebarLayout(
    sidebarPanel(
      sliderInput("threshold",
                  "Umbral:",
                  min = 0,
                  max = 1,
                  value = 0.5)
    ),
    
    mainPanel(
      plotOutput("topPlots"),
      plotOutput("bottomPlots")
    )
  )
)

# Define el servidor
server <- function(input, output) {
  
  # Matriz de referencia estática
  reference_matrix <- matrix(c(0.2, 0.5, 0.7,
                               0.1, 0.8, 0.6,
                               0.4, 0.3, 0.9), nrow = 3, byrow = TRUE)
  
  output$topPlots <- renderPlot({
    # Matriz de salida inicializada en cero
    output_matrix <- matrix(0, nrow = 3, ncol = 3)
    
    # Actualiza la matriz de salida en base al umbral
    output_matrix[reference_matrix >= input$threshold] <- "P" # Rojo (Presence)
    output_matrix[reference_matrix < input$threshold] <- "A" # Verde (Absence)
    
    # Convertir la matriz en un data frame para ggplot
    output_df <- as.data.frame(as.table(output_matrix))
    names(output_df) <- c("Row", "Col", "Value")
    
    # Convertir la matriz de referencia en un data frame para ggplot
    reference_df <- as.data.frame(as.table(reference_matrix))
    names(reference_df) <- c("Row", "Col", "RefValue")
    
    # Unir ambos data frames
    plot_df <- merge(output_df, reference_df, by = c("Row", "Col"))
    
    # Convertir las columnas de filas y columnas a factores
    plot_df$Row <- as.factor(plot_df$Row)
    plot_df$Col <- as.factor(plot_df$Col)
    
    # Definir colores y etiquetas
    color_palette <- c("P" = "red", "A" = "green")
    label_palette <- c("P" = "Presence", "A" = "Absence")
    
    # Crear el primer plot (raster plot)
    p1 <- ggplot(plot_df, aes(x = Col, y = Row, fill = Value, label = Value)) +
      geom_tile() +
      geom_text(size = 20, color = "black") +
      scale_fill_manual(values = color_palette,
                        labels = label_palette,
                        guide = guide_legend(title = "Legend")) +
      scale_y_discrete(limits = rev(levels(plot_df$Row))) +  # Invertir el eje Y
      theme_minimal() +
      labs(title = "Classified Map",
           x = "Column",
           y = "Row") +
      theme(axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            plot.title = element_text(size = 28, face = "bold", hjust = 0.5, vjust = 1),
            legend.title = element_text(size = 20, hjust = 0, vjust = 5),
            legend.text = element_text(size = 24))
    
    # Matriz adicional estática proporcionada
    static_matrix <- matrix(c("P", "A", "A",
                               "A", "P", "A",
                               "A", "A", "A"), nrow = 3, byrow = TRUE)
    
    # Convertir la matriz estática en un data frame para ggplot
    static_df <- as.data.frame(as.table(static_matrix))
    names(static_df) <- c("Row", "Col", "Value")
    
    # Convertir las columnas de filas y columnas a factores
    static_df$Row <- as.factor(static_df$Row)
    static_df$Col <- as.factor(static_df$Col)
    
    # Definir colores y etiquetas
    static_color_palette <- c("P" = "red", "A" = "green")
    static_label_palette <- c("P" = "Presence", "A" = "Absence")
    
    # Crear el segundo plot (static matrix plot)
    p2 <- ggplot(static_df, aes(x = Col, y = Row, fill = Value, label = Value)) +
      geom_tile() +
      geom_text(size = 20, color = "black") +
      scale_fill_manual(values = static_color_palette,
                        labels = static_label_palette,
                        guide = guide_legend(title = "Legend")) +
      scale_y_discrete(limits = rev(levels(static_df$Row))) +  # Invertir el eje Y
      theme_minimal() +
      labs(title = "Observed Map",
           x = "Column",
           y = "Row") +
      theme(axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            plot.title = element_text(size = 28, face = "bold", hjust = 0.5, vjust = 1),
            legend.title = element_text(size = 20, hjust = 0, vjust = 5),
            legend.text = element_text(size = 24))
    
    # Combinar los dos plots en una fila
    grid.arrange(p1, p2, nrow = 1)
  })
  
  output$bottomPlots <- renderPlot({
    # Matriz de salida inicializada en cero
    output_matrix <- matrix(0, nrow = 3, ncol = 3)
    
    # Actualiza la matriz de salida en base al umbral
    output_matrix[reference_matrix >= input$threshold] <- "P"
    output_matrix[reference_matrix < input$threshold] <- "A"
    
    # Matriz adicional estática proporcionada
    static_matrix <- matrix(c("P", "A", "A",
                              "A", "P", "A",
                              "A", "A", "A"), nrow = 3, byrow = TRUE)
    
    # Crear la nueva matriz basada en las reglas
    comparison_matrix <- matrix(NA, nrow = 3, ncol = 3)
    for (i in 1:3) {
      for (j in 1:3) {
        if (output_matrix[i, j] == "P" && static_matrix[i, j] == "P") {
          comparison_matrix[i, j] <- "Hit"
        } else if (output_matrix[i, j] == "A" && static_matrix[i, j] == "P") {
          comparison_matrix[i, j] <- "Miss"
        } else if (output_matrix[i, j] == "P" && static_matrix[i, j] == "A") {
          comparison_matrix[i, j] <- "False Alarm"
        } else if (output_matrix[i, j] == "A" && static_matrix[i, j] == "A") {
          comparison_matrix[i, j] <- "CA"
        }
      }
    }
    
    # Convertir la matriz de referencia en un data frame para ggplot
    reference_df <- as.data.frame(as.table(reference_matrix))
    names(reference_df) <- c("Row", "Col", "RefValue")
    
    # Crear el plot de la matriz de referencia
    p3 <- ggplot(reference_df, aes(x = Col, y = Row, fill = RefValue, label = round(RefValue, 2))) +
      geom_tile() +
      geom_text(size = 20, color = "black") +
      scale_fill_gradient(low = "white", high = "blue") +
      scale_y_discrete(limits = rev(levels(reference_df$Row))) +  # Invertir el eje Y
      theme_minimal() +
      labs(title = "Index Variable",
           x = "Column",
           y = "Row") +
      theme(axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            plot.title = element_text(size = 28, face = "bold", hjust = 0.5, vjust = 1),
              legend.title = element_text(size = 20, hjust = 0, vjust = 5),
              legend.text = element_text(size = 24))
    
    # Convertir la matriz de comparación en un data frame para ggplot
    comparison_df <- as.data.frame(as.table(comparison_matrix))
    names(comparison_df) <- c("Row", "Col", "Value")
    
    # Convertir las columnas de filas y columnas a factores
    comparison_df$Row <- as.factor(comparison_df$Row)
    comparison_df$Col <- as.factor(comparison_df$Col)
    
    # Crear el plot de comparación
    p4 <- ggplot(comparison_df, aes(x = Col, y = Row, fill = Value, label = Value)) +
      geom_tile() +
      geom_text(size = 6, color = "black") +
      scale_fill_manual(values = c("Hit" = "green", "Miss" = "red", 
                                   "False Alarm" = "yellow", "CA" = "lightblue")) +
      scale_y_discrete(limits = rev(levels(comparison_df$Row))) +  # Invertir el eje Y
      theme_minimal() +
      labs(title = "Confusion Matrix",
           x = "Column",
           y = "Row") +
      theme(axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            plot.title = element_text(size = 28, face = "bold", hjust = 0.5, vjust = 1),
            legend.title = element_text(size = 20, hjust = 0, vjust = 5),
            legend.text = element_text(size = 24))
    
    # Combinar los dos plots en una fila
    grid.arrange(p3, p4, nrow = 1)
  })
}

# Corre la aplicación
shinyApp(ui = ui, server = server)
