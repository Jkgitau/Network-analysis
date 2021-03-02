### This code details how to prepare the metadata file in readiness for analysis in R later on.

# Load SRA metadata from Richard et al 2016 study
metadata <- read.csv(here::here("SraRunTable.csv"),header = TRUE)

# obtain sample metadata to be used later in analysis in R.
matches <- c("Run","serologic_response_status")
samples.metadata <- metadata[grepl(paste(matches, collapse="|"), names(metadata))]


# create a grouping factor that will place each sample in the one of three vaccination states i.e unvaccinated(UV), 
# Vaccinated not protected(VNP) and vaccinated protected(VP)
status <- factor(c("UV","UV","UV","UV","UV","VNP","VNP","VNP","VNP","VNP","VNP","VNP","VNP","VNP","VNP","VNP","VNP",
                   "VNP","VNP","VNP","VNP","VNP","VNP","VNP","VNP","VNP","VNP","VNP","VNP","VP","VNP","VP","VNP",
                   "VNP","VP","VP","VP","VP","VP","VP","VP","VP","VP","VP","VP"))
                 
# append factor to samples.metadata to group samples
samples.metadata["Status"] <- status


# save the metadata to an R object
saveRDS(samples.metadata, here::here("samples.metadata.RDS"))
