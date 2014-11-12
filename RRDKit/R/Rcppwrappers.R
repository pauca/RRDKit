# This file is part of RRDKit.
# 
# RRDKit is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# Foobar is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with Foobar.  If not, see <http://www.gnu.org/licenses/>.

library(Rcpp) 

showmol<-function( ptr , open = T ){
  svg <- mol2svg(ptr)
  fileName <- tempfile(pattern = "", tmpdir = tempdir(), fileext = ".svg")
  file <- file(fileName,"w")
  cat(svg,file=file)
  close(file)
  if(open){ browseURL(paste("file:///",fileName ,sep=""))}
  return(fileName)
}

showmols<-function( ptr , open = T){
  #svgs <- molSupplierApply(ptr,showmol,open=F)
  svgs <- sapply(unlist(ptr),showmol,open=F) 
  head <- ' <!DOCTYPE html><html>
  <head>
  <style> 
    .molimg {
      height:200px
    }
    th {
      text-align:left;
    }
  </style>
  </head>  
  <body><table><tr><th>Index</th><th>Mol</th></td>'
  tail <- '</table></body></html>'

  fileName <- tempfile(pattern = "", tmpdir = tempdir(), fileext = ".html")
  file <- file(fileName,"w")
writeLines(head,file)
d <- sapply(1:length(svgs),function(i)writeLines(paste(
  '<tr><td>',i,'</td><td><img class=molimg src="',svgs[[i]], '" alt="tanding on oval"></td></tr>',sep=""),
  file))
writeLines(tail,file)
close(file)
  
 if(open){ browseURL(paste("file:///",fileName ,sep=""))}
  return(fileName)
}

showmols.grid<-function( ptr , open = T){
  svgs <- sapply(unlist(ptr),showmol,open=F)
  head <- ' <!DOCTYPE html><html>
<head>
<style>
.molimg {
height:200px
}
.molgrid {
-webkit-column-count: 5; /* Chrome, Safari, Opera */
-moz-column-count: 5; /* Firefox */
column-count: 5;
}
.molbox{
display:table;
}
</style>
</head>
<body><div class="molgrid">'
  tail <- '</div></body></html>'
  fileName <- tempfile(pattern = "", tmpdir = tempdir(), fileext = ".html")
  file <- file(fileName,"w")
  writeLines(head,file)
  d <- sapply(1:length(svgs),function(i)writeLines(paste(
    '<div class="molbox"><h3>',i,'</h3> <img class="molimg" src="',svgs[[i]], '" alt="tanding on oval"></div>',sep=""),
    file))
  writeLines(tail,file)
  close(file)
  if(open){ browseURL(paste("file:///",fileName ,sep=""))}
  return(fileName)
}


# showmols.shiny<-function( ptr , data = data.frame()){
#   require(shiny)
#    svgs <- sapply(unlist(ptr),mol2svg) 
#                     
#    svgs <- sapply( svgs, function(s) paste( "<div class='molma'>", 
#                                             gsub("svg:","",s,fixed=T),
#                                             "</div>"))
#   #svgs <- sapply(unlist(ptr),function(s) paste( "<img src='",showmol(s,F),"' />",sep="") )
#   server <- function(input, output) {
#     output$genericTable <-  renderTable({
#       if(nrow(data)==length(svgs)){
#        # cbind(data.frame(id=1:length(svgs), img=svgs),data)               
#       }else{
#         data.frame(id=1:length(svgs), img=svgs)
#       }
#     } , sanitize.text.function = function(x) x)
#   }  
#   ui <- shinyUI(fluidPage(
#     tags$head(
#       tags$style(HTML("
#        .molma > svg{
#           transform: scale(0.2, 0.2);
#        }
#        #genericTable > table > tbody > tr > :first-child{
#           display:none;
#        }
#     "))
#     ),                      
#     mainPanel(
#     tableOutput("genericTable")
#   )))
#   message("Press ESC to continue.")
#   shinyApp(ui = ui, server = server)
# }



molGetProps <- function( m ){
  if(is.list(m)){ warning("Invalid input")}
  l <- list()
  for( p in p_molGetPropList(m)){
     l[[p]]<- tryCatch( p_molGetProp(m,p), error=function(e){return (NA)})
  }
  return(l)
}

molsGetProps <- function( ms ){
  if(!is.list(ms)){ warning("Invalid input")}
  lprops<- lapply(ms,molGetProps)
  n <- sort(unique(unlist(lapply(lprops,function(l)names(l)))))
  m <- lapply(lprops,function(x) { u <- unlist(x)[n];names(u)<-n;return(u)})
  as.data.frame(do.call(rbind,m))
}

molSetProp <- function( m  ,key , value){
  p_molSetProp(m, as.character(value), key)
}

molsSetProps <- function( ms, key, v  ){
  for(i in 1:length(ms)){
    molSetProp(ms[[i]],key,v[i])
  }
}

read.sdf<-function(file){
  ms <- molSupplier(file)
  obj<- unlist(molSupplierApply(ms,function(m)return(m)))
  obj
}

write.sdf <- function(file,mols){
  p_writeSdf(file,mols)
}

