#!/opt/ohpc/pub/libs/gnu/R_base/3.3.3/bin/Rscript

# Plots the trinucleotide context (TNC) of DNA artifactual damage detected by Mutect2.
# Damage is either true negative (correctly filtered as orientation) or false positive (PASS).

library(ggplot2)
args <- commandArgs(trailingOnly = TRUE)

# Make sure caller has supplied a TNC file.
if (length(args)==0) {
  stop("Usage. Include file containing TNCs of each SBS record in caller VCF as first arguement", call.=FALSE)
}

# Get the mutational matrix for the SNVs within this allele frequence range.

# Give it one of every mutational context to start off.
# This will prevent NA's being introduced down the line if a freq range that
# does not contain any mutations of a given context is selected.
# We will correct for every context having an additional variant below before returning.

ctx = c("A[C>A]A", "A[C>A]C", "A[C>A]G", "A[C>A]T", "C[C>A]A", "C[C>A]C", "C[C>A]G", 
        "C[C>A]T", "G[C>A]A", "G[C>A]C", "G[C>A]G", "G[C>A]T", "T[C>A]A", "T[C>A]C", 
        "T[C>A]G", "T[C>A]T", "A[C>G]A", "A[C>G]C", "A[C>G]G", "A[C>G]T", "C[C>G]A", 
        "C[C>G]C", "C[C>G]G", "C[C>G]T", "G[C>G]A", "G[C>G]C", "G[C>G]G", "G[C>G]T", 
        "T[C>G]A", "T[C>G]C", "T[C>G]G", "T[C>G]T", "A[C>T]A", "A[C>T]C", "A[C>T]G", 
        "A[C>T]T", "C[C>T]A", "C[C>T]C", "C[C>T]G", "C[C>T]T", "G[C>T]A", "G[C>T]C", 
        "G[C>T]G", "G[C>T]T", "T[C>T]A", "T[C>T]C", "T[C>T]G", "T[C>T]T", "A[T>A]A", 
        "A[T>A]C", "A[T>A]G", "A[T>A]T", "C[T>A]A", "C[T>A]C", "C[T>A]G", "C[T>A]T", 
        "G[T>A]A", "G[T>A]C", "G[T>A]G", "G[T>A]T", "T[T>A]A", "T[T>A]C", "T[T>A]G", 
        "T[T>A]T", "A[T>C]A", "A[T>C]C", "A[T>C]G", "A[T>C]T", "C[T>C]A", "C[T>C]C", 
        "C[T>C]G", "C[T>C]T", "G[T>C]A", "G[T>C]C", "G[T>C]G", "G[T>C]T", "T[T>C]A", 
        "T[T>C]C", "T[T>C]G", "T[T>C]T", "A[T>G]A", "A[T>G]C", "A[T>G]G", "A[T>G]T", 
        "C[T>G]A", "C[T>G]C", "C[T>G]G", "C[T>G]T", "G[T>G]A", "G[T>G]C", "G[T>G]G", 
        "G[T>G]T", "T[T>G]A", "T[T>G]C", "T[T>G]G", "T[T>G]T")


# Get a summary breakdown.
# Read in the TNC of every variant under analysis. 
tmp = as.character(unlist(read.table(args[1], header = F))) 

# Orient contexts wrt transitions from C or T  
for(i in tmp)
{
  splits <- strsplit(i,"")[[1]]
  if(!(splits[3]=="C" || splits[3]=="T"))
  {
    tmp=chartr('ATGC', 'TACG', splits[1])
    splits[1]=chartr('ATGC', 'TACG', splits[7])
    splits[7]=tmp
    splits[3]=chartr('ATGC', 'TACG', splits[3])
    splits[5]=chartr('ATGC', 'TACG', splits[5])
    i <- paste(splits, collapse = "")
  }
}
ctx = c(ctx, tmp)

tmp=summary(as.factor(ctx))

# Mutational Order
# The order in which they will be displayed in a plot.
SIGNATURE_ORDER = c("A[C>A]A", "A[C>A]C", "A[C>A]G", "A[C>A]T", "C[C>A]A", "C[C>A]C", "C[C>A]G", 
                    "C[C>A]T", "G[C>A]A", "G[C>A]C", "G[C>A]G", "G[C>A]T", "T[C>A]A", "T[C>A]C", 
                    "T[C>A]G", "T[C>A]T", "A[C>G]A", "A[C>G]C", "A[C>G]G", "A[C>G]T", "C[C>G]A", 
                    "C[C>G]C", "C[C>G]G", "C[C>G]T", "G[C>G]A", "G[C>G]C", "G[C>G]G", "G[C>G]T", 
                    "T[C>G]A", "T[C>G]C", "T[C>G]G", "T[C>G]T", "A[C>T]A", "A[C>T]C", "A[C>T]G", 
                    "A[C>T]T", "C[C>T]A", "C[C>T]C", "C[C>T]G", "C[C>T]T", "G[C>T]A", "G[C>T]C", 
                    "G[C>T]G", "G[C>T]T", "T[C>T]A", "T[C>T]C", "T[C>T]G", "T[C>T]T", "A[T>A]A", 
                    "A[T>A]C", "A[T>A]G", "A[T>A]T", "C[T>A]A", "C[T>A]C", "C[T>A]G", "C[T>A]T", 
                    "G[T>A]A", "G[T>A]C", "G[T>A]G", "G[T>A]T", "T[T>A]A", "T[T>A]C", "T[T>A]G", 
                    "T[T>A]T", "A[T>C]A", "A[T>C]C", "A[T>C]G", "A[T>C]T", "C[T>C]A", "C[T>C]C", 
                    "C[T>C]G", "C[T>C]T", "G[T>C]A", "G[T>C]C", "G[T>C]G", "G[T>C]T", "T[T>C]A", 
                    "T[T>C]C", "T[T>C]G", "T[T>C]T", "A[T>G]A", "A[T>G]C", "A[T>G]G", "A[T>G]T", 
                    "C[T>G]A", "C[T>G]C", "C[T>G]G", "C[T>G]T", "G[T>G]A", "G[T>G]C", "G[T>G]G", 
                    "G[T>G]T", "T[T>G]A", "T[T>G]C", "T[T>G]G", "T[T>G]T")


# Match the order of the mutation types to order they will be displayed on plot.
new_order = match(SIGNATURE_ORDER, names(tmp))

# Reorder signatures vector.
tmp = tmp[new_order]
names(tmp) <- SIGNATURE_ORDER

# Correct for the extra mutation of each context we added above.
tmp = tmp-1

# Create the mutation matrix.
new_mat = matrix(tmp, 96, byrow=TRUE)
rownames(new_mat) <- names(tmp)
colnames(new_mat) <- "value"

# Sum new_max to see how many SNVs we got..(make sure to do this before normalising!!)
total_snvs = sum(new_mat)

# Normalize.
#if(total_snvs > 0)
#  new_mat = new_mat/sum(new_mat)

# Format the data and scale for the plot.
max_new_mat = max(new_mat)

# Work out y-axis scaling (roughly)
if(max_new_mat < 0.1) {
  ymax = 0.1
} else if (max_new_mat < 0.2) {
  ymax = 0.2
} else {
  ymax = 0.3
}

colors = c(
  "#02bcedff", "#000000ff", "#e22a27ff",
  "#ccc8c9ff", "#a0cf62ff", "#e7c5c3ff")

context = c(rep(c(
  "ACA", "ACC", "ACG", "ACT",
  "CCA", "CCC", "CCG", "CCT",
  "GCA", "GCC", "GCG", "GCT",
  "TCA", "TCC", "TCG", "TCT"), 3),
  rep(c(
    "ATA", "ATC", "ATG", "ATT",
    "CTA", "CTC", "CTG", "CTT",
    "GTA", "GTC", "GTG", "GTT",
    "TTA", "TTC", "TTG", "TTT"), 3))

substitution = rep(c('C>A','C>G','C>T','T>A','T>C','T>G'), each = 16)
substring(context, 2, 2) = "."
df = data.frame(substitution = substitution, context = context)
rownames(new_mat) = NULL
df1 = data.frame(substitution = substitution, context = context)
rownames(new_mat) = NULL
df2 = cbind(df1, as.data.frame(new_mat))

# Plot it.
p = ggplot(data = df2, aes(x = context, y = value, 
                           fill = substitution, width = 1)) + geom_bar(stat = "identity", 
                                                                       colour = "black", size = 0.2) +
  ggtitle(paste0(total_snvs, " SNVs.")) +
  scale_fill_manual(values = colors) + 
  facet_grid(. ~ substitution) + ylab("number of variants") + scale_y_continuous() +
  theme_bw() + 
  theme(strip.placement = "outside", strip.text = element_text(size = 12, vjust = 0, color = "black"),strip.background =element_rect(fill="#f8f8f8ff"),
        axis.title.y = element_text(size = 12, vjust = 1), 
        axis.text.y = element_text(size = 8), axis.title.x = element_text(size = 12), 
        axis.text.x = element_text(size = 5, angle = 90, 
                                   vjust = 0.4), strip.text.x = element_text(size = 9), 
        strip.text.y = element_text(size = 9), panel.grid.major.x = element_blank(), 
        panel.spacing.x = unit(0, "lines")) +
 theme(panel.background = element_blank())

pdf(file = paste0(args[1],".pdf"), width = 7, height = 5)
p
dev.off()
