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
#' hello(fname="Matthew",lname="Wipperman")
#' hello(fname="Your",lname="Name")

hello <- function(firstname, lastname) {
  cat(paste("Hello",firstname,lastname,"!"))
  system(paste("say Hello",firstname,lastname,"!"))
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

#' library("phyloseq")
#' library("data.table")
#' library("ggplot2")
#'
#' Functions used for quick plotting (from Phyloseq developers)
#'
#' @export fast_melt
#' @examples fast_melt
fast_melt = function(physeq){
  # supports "naked" otu_table as `physeq` input.
  otutab = as(otu_table(physeq), "matrix")
  if(!taxa_are_rows(physeq)){otutab <- t(otutab)}
  otudt = data.table(otutab, keep.rownames = TRUE)
  setnames(otudt, "rn", "taxaID")
  # Enforce character taxaID key
  otudt[, taxaIDchar := as.character(taxaID)]
  otudt[, taxaID := NULL]
  setnames(otudt, "taxaIDchar", "taxaID")
  # Melt count table
  mdt = melt.data.table(otudt,
                        id.vars = "taxaID",
                        variable.name = "SampleID",
                        value.name = "count")
  # Remove zeroes, NAs
  mdt <- mdt[count > 0][!is.na(count)]
  # Calculate relative abundance
  mdt[, RelativeAbundance := count / sum(count), by = SampleID]
  if(!is.null(tax_table(physeq, errorIfNULL = FALSE))){
    # If there is a tax_table, join with it. Otherwise, skip this join.
    taxdt = data.table(as(tax_table(physeq, errorIfNULL = TRUE), "matrix"), keep.rownames = TRUE)
    setnames(taxdt, "rn", "taxaID")
    # Enforce character taxaID key
    taxdt[, taxaIDchar := as.character(taxaID)]
    taxdt[, taxaID := NULL]
    setnames(taxdt, "taxaIDchar", "taxaID")
    # Join with tax table
    setkey(taxdt, "taxaID")
    setkey(mdt, "taxaID")
    mdt <- taxdt[mdt]
  }
  return(mdt)
}

#'
#' Functions used for quick plotting (from Phyloseq developers)
#'
#' @export summarize_taxa
#' @examples summarize_taxa
summarize_taxa = function(physeq, Rank, GroupBy = NULL){
  Rank <- Rank[1]
  if(!Rank %in% rank_names(physeq)){
    message("The argument to `Rank` was:\n", Rank,
            "\nBut it was not found among taxonomic ranks:\n",
            paste0(rank_names(physeq), collapse = ", "), "\n",
            "Please check the list shown above and try again.")
  }
  if(!is.null(GroupBy)){
    GroupBy <- GroupBy[1]
    if(!GroupBy %in% sample_variables(physeq)){
      message("The argument to `GroupBy` was:\n", GroupBy,
              "\nBut it was not found among sample variables:\n",
              paste0(sample_variables(physeq), collapse = ", "), "\n",
              "Please check the list shown above and try again.")
    }
  }
  # Start with fast melt
  mdt = fast_melt(physeq)
  if(!is.null(GroupBy)){
    # Add the variable indicated in `GroupBy`, if provided.
    sdt = data.table(SampleID = sample_names(physeq),
                     var1 = get_variable(physeq, GroupBy))
    setnames(sdt, "var1", GroupBy)
    # Join
    setkey(sdt, SampleID)
    setkey(mdt, SampleID)
    mdt <- sdt[mdt]
  }
  # Summarize
  summarydt = mdt[, list(meanRA = mean(RelativeAbundance),
                         sdRA = sd(RelativeAbundance),
                         minRA = min(RelativeAbundance),
                         maxRA = max(RelativeAbundance)),
                  by = c(Rank, GroupBy)]
  return(summarydt)
}

#'
#'
#' Functions used for quick plotting (from Phyloseq developers)
#'
#' @export plot_taxa_summary
#' @examples plot_taxa_summary
plot_taxa_summary = function(physeq, Rank, GroupBy = NULL){
  # Get taxa summary table
  dt1 = summarize_taxa(physeq, Rank = Rank, GroupBy = GroupBy)
  # Set factor appropriately for plotting
  RankCol = which(colnames(dt1) == Rank)
  setorder(dt1, -meanRA)
  dt1[, RankFac := factor(dt1[[Rank]],
                          levels = rev(dt1[[Rank]]))]
  dt1[, ebarMax := max(c(0, min(meanRA + sdRA))), by = eval(Rank)]
  dt1[, ebarMin := max(c(0, min(meanRA - sdRA))), by = eval(Rank)]
  # Set zeroes to one-tenth the smallest value
  ebarMinFloor = dt1[(ebarMin > 0), min(ebarMin)]
  ebarMinFloor <- ebarMinFloor / 10
  dt1[(ebarMin == 0), ebarMin := ebarMinFloor]

  pRank = ggplot(dt1, aes(x = meanRA, y = RankFac)) +
    scale_x_log10() +
    xlab("Mean Relative Abundance") +
    ylab(Rank) +
    theme_bw()
  if(!is.null(GroupBy)){
    # pRank <- pRank + facet_wrap(facets = as.formula(paste("~", GroupBy)))
    pRank <- pRank + geom_point(mapping = aes_string(colour = GroupBy),
                                size = 5)  + geom_errorbarh(aes(xmax = ebarMax,
                                                                xmin = ebarMin))
  } else {
    # Don't include error bars for faceted version
    pRank <- pRank + geom_errorbarh(aes(xmax = ebarMax,
                                        xmin = ebarMin))
  }
  return(pRank)
}
