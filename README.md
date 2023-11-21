# LECdata
repository of data for LEC plant-microbial study 
Included are files to make two phyloseq objects in R, 1 for fungal ASVS and 1 for bacterial AVS
Fungal phyloseq files include:
  Table= LECFungi.qza
  Taxonomy = merged-taxaFungiNew.qza
  Metadata = LEC-metadata4.txt
Bacterial phyloseq files include: 
  Table = table_LEC16S.qza
  Taxonomy = taxonomy_LEC16S.qza
  Rooted Tree (optional) = rooted-tree_LEC16S.qza
  Metadata = LEC-metadata4.txt
There are three r script files, one for basic microbial analyses for fungi and one for basic microbial analyses for bacteria, as well as one for more detailed analyses which can be applied to both
  Fungal R script = LECFungi.R
  Bacterial R script = LECNEW16s.R
  Other analyses = LEC_detailedanalyses.R
Vegetative data are also included in the following files: 
  Raw greenhouse data = Greenhouse_dataall.csv
  leaf number data = Numberleaf_data.csv
  Soil characteristics = SoilVariables.csv
  R script to analyze greenhouse growth data and other veg data = Greenhousedata.R
