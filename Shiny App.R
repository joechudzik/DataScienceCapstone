library(mapsapi)
library(leaflet)
library(shiny)

# Bring in redefined hclust algorithm
source('~/Documents/GitHub/DataScienceCapstone/hclust ideas/rclust.R')

# Bring in dissimilarity matrices (time and distance) to be used in rhclust
pdis <- read.csv('https://raw.githubusercontent.com/joechudzik/DataScienceCapstone/master/Data/pdist.csv')
row.names(pdis) <- pdis[,1]; pdis <- pdis[,-1]
ptim <- read.csv('https://raw.githubusercontent.com/joechudzik/DataScienceCapstone/master/Data/ptime.csv')
row.names(ptim) <- ptim$X; ptim <- ptim[,-1]

# Main data file
data <- read.csv('~/Documents/GitHub/DataScienceCapstone/Data/simulatedData.csv')
data <- data[1:10,]

# Vectors for origins and destinations.
origin_vector <- paste0(data$HSE_NBR_home,' ', data$STREET_home, ' ', data$STTYPE_home, ', Milwaukee, WI ', data$ZIP_CODE_home)
destination_vector <- paste0(data$HSE_NBR,' ', data$STREET, ' ', data$STTYPE, ', Milwaukee, WI ', data$ZIP_CODE)

doc=mp_directions(
  origin = origin_vector[1],
  destination = destination_vector[1],
  alternatives=FALSE,
  key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
)

routes <- mp_get_routes(doc)

for(x in 2:length(origin_vector)){
  doc <- mp_directions(
    origin = origin_vector[x],
    destination = destination_vector[x],
    alternatives = FALSE,
    key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
  )
  route <- mp_get_routes(doc)
  
  routes <- rbind(route, routes)
}

r_colors <- rgb(t(col2rgb(colors()) / 255))
names(r_colors) <- colors()

ui <- fluidPage(
  leafletOutput("mymap"),
  p(),
    titlePanel("Actions"),
    
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar to demonstrate various slider options ----
    sidebarPanel(
      
      # Input: Simple integer interval ----
      sliderInput("timeImportanceWeight", "Time Factor:",
                  min = 0, max = 1,
                  value = 0.5),
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Table summarizing the values entered ----
      tableOutput("values"),
      plotOutput('optimalDendrogram'),
      textOutput('optimalCompletePath')
    )
  )
)
pal <- colorFactor(palette='Dark2', domain=rownames(routes))
server <- function(input, output, session) {
  
  points <- eventReactive(input$recalc, {
    cbind(rnorm(40) * 2 + 13, rnorm(40) + 48)
  }, ignoreNULL = FALSE)
  
  output$mymap <- renderLeaflet({
    leaflet() %>%
      addProviderTiles("CartoDB.Positron") %>%
      addPolylines(data = routes[1], opacity = 1, weight = 7, color = ~pal(alternative_id))
  })
  
  sliderValues <- reactive({
    
    data.frame(
      Name = c("Time Factor"),
      Value = as.character(c(input$timeImportanceWeight)),
      stringsAsFactors = FALSE)
    
  })
  
  output$optimalDendrogram <- renderPlot({
    optimalCluster <- rhclust(pdis, ptim, input$timeImportanceWeight)
    optimalPlot <- plot(optimalCluster); optimalPlot <- rect.hclust(optimalCluster, k=3, border=2:6)
    optimalPlot
  })
  
  output$optimalCompletePath <- renderText({
    optimalCluster <- rhclust(pdis, ptim, input$timeImportanceWeight)
    optimalPath <- optimalCluster$merge.route
    paste0(optimalPath[[length(optimalPath)]])
  })
  
  # Show the values in an HTML table ----
  output$values <- renderTable({
    sliderValues()
  })
}

shinyApp(ui, server)
