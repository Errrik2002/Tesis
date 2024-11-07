install.packages("shiny")
library(shiny)
library(bslib)


#Define UI
ui <- page_sidebar(
  title = "tittle panel",
  sidebar = sidebar("sidebar"),
  "main content",
  sliderInput(
    "range",
    label = "Set value",
    min = 0,
    max= 100,
    value= c(0,100)
  )
)

#degine server logic
server <- function(input, output){
  
}

#Run the APP
shinyApp(ui = ui, server = server)



