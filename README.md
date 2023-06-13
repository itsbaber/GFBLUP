# GFBLUP
Genomic Features Unbiased Linear Prediction (GFBLUP)
GFBLUP is a strategy to incorporate biological information in genomic prediciton models. SNP markers can be categorized based on a prior biological information. The information can originate from different biological sources like Gene Ontology information, known causal variants, omics data, and every piece of information that can be linked to a set of SNPs. 
The categorization requires three things including Chromosome number, start and end positions, and of course SNPs (markers) themselve.
In GFBLUP, we generally categorize SNPs into two categories, selected and remaining. However, even categories are possible based on biological importance or anyother supported information. 
GFBLUP considers selected and remaining categories to be independent from each other, i.e., they don't share SNPs with each other. However, if we are using multiple pieces of biological information, it means we can perform GFBLUP for each of the biological feature, and in that cases SNPs can be shared between different features; however, selected and remaining SNPs categories remain independent from each other.

The R based algorithm to perform SNP categorization can be found in the file "SNP_cats.R"
