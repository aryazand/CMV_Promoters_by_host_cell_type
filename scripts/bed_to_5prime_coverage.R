make_fiveprimeends_grange <- function(input_bed) {
  
  fiveprime <- input_bed |> GenomicRanges::resize(width = 1, fix = "start")
  fiveprime_unique <- fiveprime |> BiocGenerics::unique()
  fiveprime_unique_counts <- GenomicRanges::countOverlaps(fiveprime_unique, fiveprime)
  BiocGenerics::score(fiveprime_unique) <- fiveprime_unique_counts
  
  return(fiveprime_unique)
}

x = make_fiveprimeends_grange(snakemake@input[[1]]) 
rtracklayer::export.bed(x, snakemake@output[[1]])