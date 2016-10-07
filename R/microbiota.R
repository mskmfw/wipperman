#' @import plyr
#' @import dplyr
#' @import phyloseq
#' @import yingtools2
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

#' format.lefse.table
#'
#' See Curtis Huttenhower's LeFSe
#'
#' @export
#' @examples format.lefse.table takes a tax.file from blastn output (either 16SMicrobial or refseq_rna) and formats lefse table with class "in quotes
lefse.format <- function(phyloseq,class) {
  t <- fread(tax.file,colClasses=c("sallgi"="character","staxids"="character")) %>% tbl_df() %>%
    mutate(taxonomy=gsub("\\[(superkingdom|phylum|class|order|family|genus|species)\\]","",taxonomy),
           staxid=as.numeric(sapply(strsplit(staxids,split=";"),first)),
           otu=paste0(qseqid,";"),
           otu.number=as.numeric(str_extract(otu,"(?<=OTU_)[0-9]+"))) %>%
    separate(taxonomy,into=c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),sep="\\|",remove=FALSE) %>%
    group_by(otu) %>% arrange(evalue,staxid) %>% filter(!duplicated(taxonomy)) %>%
    mutate(n.ties=sum(dense_rank(evalue)==1),blast.data=paste0(Species," (eval=",evalue,",pid=",pident,")",collapse=";")) %>%
    filter(row_number()==1) %>% ungroup() %>% arrange(otu.number) %>%
    dplyr::select(otu,evalue,pident,Kingdom,Phylum,Class,Order,Phylum,Class,Order,Family,Genus,Species,n.ties,blast.data,everything())
  dat <- t[,c("otu","taxonomy")] #two columns with otu and taxonomy info separated by | symbol for lefse
  otus <- otu_table(phyloseq)
  x <- as.data.frame(otus)
  x <- cbind(otu = rownames(x), x)
  y <- as.data.frame(t(sample_data(phyloseq)))
  var <- y[class,]
  names <- rownames(var)
  rownames(var) <- NULL
  var <- cbind(names,var)
  colnames(var)[1] <- "taxonomy"
  #write.table(y, file="info for lefse.txt", sep="\t")
  ttt <- left_join(x,dat)
  ttt$otu <- NULL
  col_idx <- grep("taxonomy", names(ttt))
  ttt <- ttt[, c(col_idx, (1:ncol(ttt))[-col_idx])]
  ttx <- rbind(ttt,var)
  row_idx <- grep(class,ttx[,1])
  ttxx <- ttx[c(row_idx,(1:nrow(ttx))[-row_idx]),]
  write.table(ttxx, file="otus_lefse.txt", sep="\t",row.names = F)
}

#' read.otu.melt.phyloseq
#'
#' Should load with Ying's Tools (ying14/yingtools) but doesn't for some reason
#'
#' @export
#' @examples read.otu.melt.phyloseq makes a nice table
read.otu.melt.phyloseq <- function(phy,filter.zero=TRUE,sample_data=TRUE) {
  #phy0=phy;phy=subset_taxa(phy0,taxa_names(phy0) %in% head(taxa_names(phy0),10))
  otu0 <- otu_table(phy) %>% as.matrix() %>% melt(varnames=c("otu","sample"),value.name="numseqs") %>% as_data_frame() %>% mutate(otu=as.character(otu),sample=as.character(sample))
  tax0 <- get.tax(phy)
  tax0.match <- select(tax0,-otu)[match(otu0$otu,tax0$otu),]
  otu <- cbind(otu0,tax0.match) %>%
    group_by(sample) %>% mutate(pctseqs=prop.table(numseqs)) %>% ungroup() %>% tbl_df()
  if (filter.zero) {
    otu <- otu %>% filter(numseqs>0)
  }
  #add sample data
  if (sample_data) {
    samp0 <- get.samp(phy)
    otu <- otu %>% left_join(samp0,by="sample")
  }
  return(otu)
}

