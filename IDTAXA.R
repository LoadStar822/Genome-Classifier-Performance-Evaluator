library(DECIPHER)
library(pryr)

dir_path <- "/dev/disk5/zqtqz/datasets/zong"

result_dir <- "/dev/disk5/zqtqz/project/zongshu/result/idtaxa_results"

fasta_files <- list.files(path = dir_path, pattern = "\\.fasta$")

total_time <- 0
total_mem <- 0

for (file in fasta_files) {
  full_path <- file.path(dir_path, file)

  mem_before <- mem_used()

  time_taken <- system.time({
    seqs <- readDNAStringSet(full_path)
    seqs <- RemoveGaps(seqs)

    load("<<REPLACE WITH PATH TO RData file>>")

    ids <- IdTaxa(seqs, trainingSet, strand="both", threshold=60, processors=NULL)
  })

  mem_after <- mem_used()

  total_time <- total_time + time_taken[3]
  total_mem <- total_mem + (mem_after - mem_before)

  print(paste("Time taken for file", file, ":", time_taken))

  print(paste("Memory used for file", file, ":", mem_after - mem_before))

  result_file <- gsub(".fasta", ".rds", file)
  result_path <- file.path(result_dir, result_file)

  saveRDS(ids, file = result_path)

  print(ids)
  plot(ids)

  write(paste("Total time:", total_time, "Total memory:", total_mem), file = file.path(result_dir, "results.txt"), append = TRUE)
}
