#' Print hello
#'
#' This function prints hello
#'
#' @param fname First name
#' @param lname Last name
#' @export
#' @examples
#' hello(fname="Tomas",lname="Greif")
#' hello(fname="Your",lname="Name")

hello <- function(fname, lname) {
  cat(paste("Hello",fname,lname,"!"))
  system(paste("say Hello",fname,lname,"!"))
}

#' plot.lefse.results
#'
#' This function plots the results from Curtis Huttenhower's LeFSe output
#'
#' @param
#' @param
#' @export
#' @examples plot.lefse.results() #must have a file called lefse.res in working directory
#' plot.lefse.results()
plot.lefse.results <- function(){
  #reads a file called lefse.res, which is output from the system running LeFSe
  lefse.res <- read.table(file="lefse.res",sep="\t")%>%
    rename(taxon = V1, log.max.pct = V2, direction = V3,
           lda = V4, p.value = V5)%>%as.data.frame() %>%na.omit(lefse.res)%>%
    mutate(taxon=gsub("\\.","|",taxon))
  lefse.res$taxon <- factor(lefse.res$taxon, levels = lefse.res$taxon[order(lefse.res$direction,lefse.res$lda)])
  p <- ggplot(lefse.res, aes(x=taxon,y=lda,fill=direction)) +
    geom_bar(stat = "identity",width=0.1,position=position_dodge(1)) + coord_flip() +
    geom_text(aes(label=p.value), vjust=0,hjust=-0.1) +
    theme(legend.position="bottom")
  return(p)
}
