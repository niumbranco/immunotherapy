# libs
library(tidyverse)
library(RColorBrewer)
library(vegan)
library(phyloseq)
library(scales)
library(cowplot)

root_dir <- "~/Repositories/nium/metasoil"
out_dir <-  "~/Repositories/nium/metasoil/results"

dir.create(out_dir, showWarnings = F)

## load phyloseq objects
its_phyobj_all = readRDS( file.path(root_dir, "RDS_files", "ITS_all_phyobj.rds"))
#its_phyobj_10k = readRDS( here("RDS_files/ITS_rarefied10k_phyobj.rds"))

D_phyobj_all = readRDS(file.path(root_dir, "RDS_files","OTU_DNA_all_phyobj.rds"))
R_phyobj_all = readRDS(file.path(root_dir, "RDS_files", "OTU_RNA_all_phyobj.rds"))

sampleinfo = read_tsv(file.path(root_dir, "processed_tables/samplemetadata.tsv"), show_col_types = FALSE)
samplelist = read_tsv(file.path(root_dir, "processed_tables/sampleIDs_info_after_rarefaction.tsv"), show_col_types = FALSE)

samples_retain = samplelist %>% 
  filter(is.na(Remove))  %>% 
  pull(sampleID) %>%  unique()

sampleinfo = sampleinfo %>% 
  filter(sampleID %in% samples_retain)

## aesthetics
### colors for land use
#TODO generalize

# land_use_names = unique(sampleinfo$LandUse)
land_use_names = c("Arable", "Grassland", "Forest","Vineyard", "Orchard" )
land_use_cols = brewer.pal(length(land_use_names), "Set1")
land_use_cols = setNames( land_use_cols, land_use_names )

#testing#phyobj = its_phyobj_all
#############

# parallelised rarefaction curve function
#https://dave-clark.github.io/post/speeding-up-rarefaction-curves-for-microbial-community-ecology/
quickRareCurve <- function (x, step = 1, sample, xlab = "Sample Size",
                            ylab = "Species", label = TRUE, col, lty, max.cores = T, nCores = 1, ...)
{
  require(parallel)
  x <- as.matrix(x)
  if (!identical(all.equal(x, round(x)), TRUE))
    stop("function accepts only integers (counts)")
  if (missing(col))
    col <- par("col")
  if (missing(lty))
    lty <- par("lty")
  tot <- rowSums(x) # calculates library sizes
  S <- specnumber(x) # calculates n species for each sample
  if (any(S <= 0)) {
    message("empty rows removed")
    x <- x[S > 0, , drop = FALSE]
    tot <- tot[S > 0]
    S <- S[S > 0]
  } # removes any empty rows
  nr <- nrow(x) # number of samples
  col <- rep(col, length.out = nr)
  lty <- rep(lty, length.out = nr)
  # parallel mclapply
  # set number of cores
  mc <- getOption("mc.cores", ifelse(max.cores, detectCores(), nCores))
  message(paste("Using ", mc, " cores"))
  out <- mclapply(seq_len(nr), mc.cores = mc, function(i) {
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i])
      n <- c(n, tot[i])
    drop(rarefy(x[i, ], n))
  })
  Nmax <- sapply(out, function(x) max(attr(x, "Subsample")))
  Smax <- sapply(out, max)
  plot(c(1, max(Nmax)), c(1, max(Smax)), xlab = xlab, ylab = ylab,
       type = "n", ...)
  if (!missing(sample)) {
    abline(v = sample)
    rare <- sapply(out, function(z) approx(x = attr(z, "Subsample"),
                                           y = z, xout = sample, rule = 1)$y)
    abline(h = rare, lwd = 0.5)
  }
  for (ln in seq_along(out)) {
    N <- attr(out[[ln]], "Subsample")
    lines(N, out[[ln]], col = col[ln], lty = lty[ln], ...)
  }
  if (label) {
    ordilabel(cbind(tot, S), labels = rownames(x), ...)
  }
  invisible(out)
}


#get data for rarefaction curves
rarefaction_plot <- function(phyobj, RAREFACTION_CUTOFF) {
  # extract tables from phyloseq object
  OTU = otu_table(phyobj) %>% 
    as.matrix() %>% 
    as.data.frame() %>% 
    rownames_to_column("OTU")
  SAMPLES = sample_data(phyobj) %>% 
    as.matrix() %>% 
    as.data.frame() %>% 
    rownames_to_column("sampleID")
  
  #tt = rarecurve(t(OTU[-1]), step=50, cex=0.5, label=T)
  tt = quickRareCurve(t(OTU[-1]), step=50, cex=0.5, label=T) # uses all cores by default
  #https://stackoverflow.com/questions/47234809/coloring-rarefaction-curve-lines-by-metadata-vegan-package-phyloseq-package
  rare <- lapply(tt, function(x){
    b <- as.data.frame(x)
    b <- data.frame(OTU = b[,1], raw.read = rownames(b))
    b$raw.read <- as.numeric(gsub("N", "",  b$raw.read))
    return(b)
  })
  names(rare) <- rownames(t(OTU[-1]))
  rare <- map_dfr(rare, function(x){
    z <- data.frame(x)
    return(z)
  }, .id = "sampleID")
  
  rareplot = rare %>%  
    left_join(SAMPLES) %>% 
    mutate(year_season= paste(Year, Season, sep="_")) %>% 
    ggplot((aes( x = raw.read, y = OTU, group = sampleID, color = LandUse ))) +
    geom_line(alpha = 0.7) +
    #scale_x_continuous(labels =  scales::scientific_format())+
    scale_color_manual(values = land_use_cols) + 
    theme_bw() +
    xlab("Number of readcounts") +
    ylab("Number of OTUs detected") + 
    geom_vline(xintercept = RAREFACTION_CUTOFF, alpha =0.8, color = "black", linetype="dashed") + 
    facet_wrap(~ year_season, ncol=1) + 
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10))
  return(rareplot)
}

its_plot = rarefaction_plot(its_phyobj_all, 10000)
RNA_plot = rarefaction_plot(R_phyobj_all, 30000)
DNA_plot = rarefaction_plot(D_phyobj_all, 30000)

#################### save plots

pdf(file.path(out_dir, "0_rarefaction_16S_DNA.pdf"), 
    width = 18/2.54, 
    height = 20/2.54)
DNA_plot
dev.off()

pdf(file.path(out_dir, "0_rarefaction_16S_RNA.pdf"), 
    width = 18/2.54, 
    height = 20/2.54)
print(RNA_plot)
dev.off()

pdf(file.path(out_dir, "0_rarefaction_ITS.pdf"), 
    width = 18/2.54, 
    height = 20/2.54)
print(its_plot)
dev.off()

# arrange the three plots in a single row
prow <- plot_grid(
  DNA_plot + theme(legend.position="none") +
    scale_x_continuous(labels = label_number(scale_cut = cut_si(""))),
  RNA_plot + theme(legend.position="none",
                   axis.title.y = element_blank()) +
    scale_x_continuous(labels = label_number(scale_cut = cut_si(""))),
  its_plot + theme(legend.position="none",
                   axis.title.y = element_blank()) +
    scale_x_continuous(labels = label_number(scale_cut = cut_si(""))),
  align = 'vh',
  labels = c("A", "B", "C"),
  hjust = -1,
  nrow = 1
)

legend_b <- get_legend(
  DNA_plot + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom")
)

# add the legend underneath the row we made earlier. Give it 10%
# of the height of one plot (via rel_heights).
finplot = plot_grid(prow, legend_b, ncol = 1, rel_heights = c(1, .1))

pdf(file.path(out_dir, "0_rarefaction_all.pdf"), 
    width = 28/2.54, 
    height = 20/2.54)
print(finplot)
dev.off()
