#trawlDiversity

---

## Gradual changes in range size accompany long-term trends in species richness

**Authors** Ryan D. Batt, James W. Morley, Rebecca L. Selden, Morgan W. Tingley, Malin L. Pinsky

This R package is designed to reproduce the statistical analysis and visualization of a substantial portion of my postdoctoral research, which is in progress.

trawlDiversity relies on the [trawlData](https://github.com/rBatt/trawlData) R package. `trawlData` contains scientific bottom trawl surveys of marine animals all around the North American coastline during the past several decades. 

Interestingly, many of these species are known to be shifting their geographic range in response to climate change (see (OceanAdapt Website)[http://oceanadapt.rutgers.edu/] and (OceanAdapt Repo)[https://github.com/mpinsky/OceanAdapt]). However, what does this mean for *groups* of interacting species?

My goal is to understand how the players in these food webs may have changed -- have species been gained, or lost? If species richness has changed, has it happened with gradual or rapid shifts in the ranges of individual species? The answers to these questions will help us understand how these systems have changed in the past, and how marine ecosystems may respond to environmental change in the future.

A centerpiece of this analysis is the use of a multispecies occupancy model (MSOM) to estimate species richness while accounting for imperfect detection of species by the sampling methods use in the surveys (this is a common problem with any biodiversity study). These corrections have never been made in marine studies before, but can crucially effect estimates of biodiversity.

The MSOM models run in Stan or JAGS. These are additional programs that won't be installed directly through R. Furthermore, the models are computationally quite intensive. For this reason, this repository not only contains final results of the analysis, but also contains intermediate results.




