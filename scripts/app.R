server <- function(input, output) {
  #----------ALL DATA TABLE----------#
  output$alldata=renderDT({
    datatable(Table1,options=list(pageLength=10,
                                  lengthMenu=c(2,5,8,10)
    ),rownames = FALSE)
  })
  
  #---------TOP AND BOTTOM TABLES--------#
  TopRankData=reactive({
    if (input$topbottomselect == "Top 5") {
      ordered=head(Table1[order(-Table1$Calories), ],n=5)
    } else {
      ordered=head(Table1[order(Table1$Calories), ],n=5)
    }
    return(ordered)
  })
  
  BotRankData=reactive({
    if (input$topbottomselect == "Top 5") {
      ordered=head(Table1[order(-Table1$AvgHR), ],n=5)
    } else {
      ordered=head(Table1[order(Table1$AvgHR), ],n=5)
    }
    return(ordered)
  })
  
  output$topbottom1=renderDT({
    datatable(TopRankData(),rownames = FALSE,
              options = list(pageLength = 5, lengthChange = FALSE, dom='t'),
              caption = 'Top 5/Bottom 5 Classes by Calories Burned')
  })
  
  output$topbottom2=renderDT({
    datatable(BotRankData(),rownames = FALSE,
              options = list(pageLength =5, lengthChange = FALSE, dom='t'),
              caption = "Top 5/Bottom 5 Classes by Avg HR")
  })
}
