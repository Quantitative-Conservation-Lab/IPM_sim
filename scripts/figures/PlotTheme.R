# CUSTOM PLOT THEME

theme_ipmsim <- function(){ 
  font <- "Helvetica"   #assign font family up front
  
  theme_minimal() %+replace%    #replace elements we want to change
    
    theme(
      
      #grid elements
      panel.grid.major = element_blank(),    #strip major gridlines
      panel.grid.minor = element_blank(),    #strip minor gridlines
      
      #since theme_minimal() already strips axis lines, 
      #we don't need to do that again
      
      #text elements
      plot.title = element_text(             #title
        family = font, color = "black",      #set font family
        size = 20,                #set font size
        face = 'bold',            #bold typeface
        hjust = 0,                #left align
        vjust = 2),               #raise slightly
      
      plot.subtitle = element_text(          #subtitle
        family = font, color = "black",           #font family
        size = 14),               #font size
      
      axis.title = element_text(             #axis titles
        family = font, color = "black",           #font family
        size = 12),               #font size
      
      axis.ticks.x.bottom = element_line(colour = "grey20"), 
      axis.ticks.y.left = element_line(colour = "grey20"), 
      
      axis.text = element_text(              #axis text
        family = font, color = "black",           #axis family
        size = 12),                #font size
      
      legend.title = element_text(             #axis titles
        family = font, color = "black",           #font family
        size = 12),               #font size
      
      legend.text = element_text(             #axis titles
        family = font, color = "black",           #font family
        size = 12),              #font size
      
      strip.text = element_text(             #axis titles
        family = font, color = "black",           #font family
        size = 12)               #font size
      
      #since the legend often requires manual tweaking 
      #based on plot content, don't define it here
    )
}
