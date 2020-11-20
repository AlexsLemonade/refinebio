# Turn warnings into errors because biocLite throws warnings instead # of error if it fails to install something.
options(warn=2)
options(repos=structure(c(CRAN="https://cloud.r-project.org")))
options(Ncpus=parallel::detectCores())

# Bioconductor packages, installed by devtools::install_url()

# Helper function that installs a list of packages based on input URL
install_with_url <- function(main_url, packages) {
  pkg_ids <- devtools::install_url(paste0(main_url, packages))
  if(any(is.na(pkg_ids))) {
    pkg_fails <- paste(packages[is.na(pkg_ids)], collapse = "; ")
    stop(paste("Failed to install package(s):", pkg_fails ))
  }
  return(pkg_ids)
}

devtools::install_version('dplyr', version='1.0.2')

bioc_url <- 'https://bioconductor.org/packages/release/bioc/src/contrib/'
bioc_pkgs <- c(
  'oligo_1.54.1.tar.gz',
  'GEOquery_2.58.0.tar.gz',
  'SCAN.UPC_2.32.0.tar.gz',
  'affy_1.66.0.tar.gz',
  'affyio_1.60.0.tar.gz',
  'AnnotationDbi_1.52.0.tar.gz',
  'zlibbioc_1.36.0.tar.gz',
  'preprocessCore_1.50.0.tar.gz',
  'genefilter_1.70.0.tar.gz',
  'sva_3.36.0.tar.gz',
  'limma_3.46.0.tar.gz'
)
install_with_url(bioc_url, bioc_pkgs)

# Invoke another R script to install BrainArray ensg packages
source("install_ensg_pkgs.R")

# Install Bioconductor platform design (pd) packages
experiment_url <- 'https://bioconductor.org/packages/release/data/experiment/src/contrib/'
pd_experiment_pkgs <- c(
  'pd.atdschip.tiling_0.16.0.tar.gz'
)
install_with_url(experiment_url, pd_experiment_pkgs)

annotation_url <- 'https://bioconductor.org/packages/release/data/annotation/src/contrib/'
pd_annotation_pkgs <- c(
  'pd.081229.hg18.promoter.medip.hx1_0.99.4.tar.gz',
  'pd.2006.07.18.hg18.refseq.promoter_1.8.1.tar.gz',
  'pd.2006.07.18.mm8.refseq.promoter_0.99.3.tar.gz',
  'pd.2006.10.31.rn34.refseq.promoter_0.99.3.tar.gz',
  'pd.ag_3.12.0.tar.gz',
  'pd.aragene.1.0.st_3.12.0.tar.gz',
  'pd.aragene.1.1.st_3.12.0.tar.gz',
  'pd.ath1.121501_3.12.0.tar.gz',
  'pd.barley1_3.12.0.tar.gz',
  'pd.bovgene.1.0.st_3.12.0.tar.gz',
  'pd.bovgene.1.1.st_3.12.0.tar.gz',
  'pd.bovine_3.12.0.tar.gz',
  'pd.bsubtilis_3.12.0.tar.gz',
  'pd.cangene.1.0.st_3.12.0.tar.gz',
  'pd.cangene.1.1.st_3.12.0.tar.gz',
  'pd.canine_3.12.0.tar.gz',
  'pd.canine.2_3.12.0.tar.gz',
  'pd.celegans_3.12.0.tar.gz',
  'pd.charm.hg18.example_0.99.4.tar.gz',
  'pd.chicken_3.12.0.tar.gz',
  'pd.chigene.1.0.st_3.12.0.tar.gz',
  'pd.chigene.1.1.st_3.12.0.tar.gz',
  'pd.chogene.2.0.st_3.12.0.tar.gz',
  'pd.chogene.2.1.st_3.12.0.tar.gz',
  'pd.citrus_3.12.0.tar.gz',
  'pd.clariom.d.human_3.14.1.tar.gz',
  'pd.clariom.s.human_3.14.1.tar.gz',
  'pd.clariom.s.human.ht_3.14.1.tar.gz',
  'pd.clariom.s.mouse_3.14.1.tar.gz',
  'pd.clariom.s.mouse.ht_3.14.1.tar.gz',
  'pd.clariom.s.rat_3.14.1.tar.gz',
  'pd.clariom.s.rat.ht_3.14.1.tar.gz',
  'pd.cotton_3.12.0.tar.gz',
  'pd.cyngene.1.0.st_3.12.0.tar.gz',
  'pd.cyngene.1.1.st_3.12.0.tar.gz',
  'pd.cyrgene.1.0.st_3.12.0.tar.gz',
  'pd.cyrgene.1.1.st_3.12.0.tar.gz',
  'pd.cytogenetics.array_3.12.0.tar.gz',
  'pd.drogene.1.0.st_3.12.0.tar.gz',
  'pd.drogene.1.1.st_3.12.0.tar.gz',
  'pd.drosgenome1_3.12.0.tar.gz',
  'pd.drosophila.2_3.12.0.tar.gz',
  'pd.e.coli.2_3.12.0.tar.gz',
  'pd.ecoli_3.12.0.tar.gz',
  'pd.ecoli.asv2_3.12.0.tar.gz',
  'pd.elegene.1.0.st_3.12.0.tar.gz',
  'pd.elegene.1.1.st_3.12.0.tar.gz',
  'pd.equgene.1.0.st_3.12.0.tar.gz',
  'pd.equgene.1.1.st_3.12.0.tar.gz',
  'pd.feinberg.hg18.me.hx1_0.99.3.tar.gz',
  'pd.feinberg.mm8.me.hx1_0.99.3.tar.gz',
  'pd.felgene.1.0.st_3.12.0.tar.gz',
  'pd.felgene.1.1.st_3.12.0.tar.gz',
  'pd.fingene.1.0.st_3.12.0.tar.gz',
  'pd.fingene.1.1.st_3.12.0.tar.gz',
  'pd.genomewidesnp.5_3.14.1.tar.gz',
  'pd.genomewidesnp.6_3.14.1.tar.gz',
  'pd.guigene.1.0.st_3.12.0.tar.gz',
  'pd.guigene.1.1.st_3.12.0.tar.gz',
  'pd.hc.g110_3.12.0.tar.gz',
  'pd.hg.focus_3.12.0.tar.gz',
  'pd.hg.u133.plus.2_3.12.0.tar.gz',
  'pd.hg.u133a_3.12.0.tar.gz',
  'pd.hg.u133a.2_3.12.0.tar.gz',
  'pd.hg.u133a.tag_3.12.0.tar.gz',
  'pd.hg.u133b_3.12.0.tar.gz',
  'pd.hg.u219_3.12.0.tar.gz',
  'pd.hg.u95a_3.12.0.tar.gz',
  'pd.hg.u95av2_3.12.0.tar.gz',
  'pd.hg.u95b_3.12.0.tar.gz',
  'pd.hg.u95c_3.12.0.tar.gz',
  'pd.hg.u95d_3.12.0.tar.gz',
  'pd.hg.u95e_3.12.0.tar.gz',
  'pd.hg18.60mer.expr_3.12.0.tar.gz',
  'pd.ht.hg.u133.plus.pm_3.12.0.tar.gz',
  'pd.ht.hg.u133a_3.12.0.tar.gz',
  'pd.ht.mg.430a_3.12.0.tar.gz',
  'pd.hta.2.0_3.12.2.tar.gz',
  'pd.hu6800_3.12.0.tar.gz',
  'pd.huex.1.0.st.v2_3.14.1.tar.gz',
  'pd.hugene.1.0.st.v1_3.14.1.tar.gz',
  'pd.hugene.1.1.st.v1_3.14.1.tar.gz',
  'pd.hugene.2.0.st_3.14.1.tar.gz',
  'pd.hugene.2.1.st_3.14.1.tar.gz',
  'pd.maize_3.12.0.tar.gz',
  'pd.mapping250k.nsp_3.12.0.tar.gz',
  'pd.mapping250k.sty_3.12.0.tar.gz',
  'pd.mapping50k.hind240_3.12.0.tar.gz',
  'pd.mapping50k.xba240_3.12.0.tar.gz',
  'pd.margene.1.0.st_3.12.0.tar.gz',
  'pd.margene.1.1.st_3.12.0.tar.gz',
  'pd.medgene.1.0.st_3.12.0.tar.gz',
  'pd.medgene.1.1.st_3.12.0.tar.gz',
  'pd.medicago_3.12.0.tar.gz',
  'pd.mg.u74a_3.12.0.tar.gz',
  'pd.mg.u74av2_3.12.0.tar.gz',
  'pd.mg.u74b_3.12.0.tar.gz',
  'pd.mg.u74bv2_3.12.0.tar.gz',
  'pd.mg.u74c_3.12.0.tar.gz',
  'pd.mg.u74cv2_3.12.0.tar.gz',
  'pd.mirna.1.0_3.12.0.tar.gz',
  'pd.mirna.2.0_3.12.0.tar.gz',
  'pd.mirna.3.0_3.12.0.tar.gz',
  'pd.mirna.3.1_3.8.1.tar.gz',
  'pd.mirna.4.0_3.12.0.tar.gz',
  'pd.moe430a_3.12.0.tar.gz',
  'pd.moe430b_3.12.0.tar.gz',
  'pd.moex.1.0.st.v1_3.14.1.tar.gz',
  'pd.mogene.1.0.st.v1_3.14.1.tar.gz',
  'pd.mogene.1.1.st.v1_3.14.1.tar.gz',
  'pd.mogene.2.0.st_3.14.1.tar.gz',
  'pd.mogene.2.1.st_3.14.1.tar.gz',
  'pd.mouse430.2_3.12.0.tar.gz',
  'pd.mouse430a.2_3.12.0.tar.gz',
  'pd.mta.1.0_3.12.0.tar.gz',
  'pd.mu11ksuba_3.12.0.tar.gz',
  'pd.mu11ksubb_3.12.0.tar.gz',
  'pd.nugo.hs1a520180_3.4.0.tar.gz',
  'pd.nugo.mm1a520177_3.4.0.tar.gz',
  'pd.ovigene.1.0.st_3.12.0.tar.gz',
  'pd.ovigene.1.1.st_3.12.0.tar.gz',
  'pd.pae.g1a_3.12.0.tar.gz',
  'pd.plasmodium.anopheles_3.12.0.tar.gz',
  'pd.poplar_3.12.0.tar.gz',
  'pd.porcine_3.12.0.tar.gz',
  'pd.porgene.1.0.st_3.12.0.tar.gz',
  'pd.porgene.1.1.st_3.12.0.tar.gz',
  'pd.rabgene.1.0.st_3.12.0.tar.gz',
  'pd.rabgene.1.1.st_3.12.0.tar.gz',
  'pd.rae230a_3.12.0.tar.gz',
  'pd.rae230b_3.12.0.tar.gz',
  'pd.raex.1.0.st.v1_3.14.1.tar.gz',
  'pd.ragene.1.0.st.v1_3.14.1.tar.gz',
  'pd.ragene.1.1.st.v1_3.14.1.tar.gz',
  'pd.ragene.2.0.st_3.14.1.tar.gz',
  'pd.ragene.2.1.st_3.14.1.tar.gz',
  'pd.rat230.2_3.12.0.tar.gz',
  'pd.rcngene.1.0.st_3.12.0.tar.gz',
  'pd.rcngene.1.1.st_3.12.0.tar.gz',
  'pd.rg.u34a_3.12.0.tar.gz',
  'pd.rg.u34b_3.12.0.tar.gz',
  'pd.rg.u34c_3.12.0.tar.gz',
  'pd.rhegene.1.0.st_3.12.0.tar.gz',
  'pd.rhegene.1.1.st_3.12.0.tar.gz',
  'pd.rhesus_3.12.0.tar.gz',
  'pd.rice_3.12.0.tar.gz',
  'pd.rjpgene.1.0.st_3.12.0.tar.gz',
  'pd.rjpgene.1.1.st_3.12.0.tar.gz',
  'pd.rn.u34_3.12.0.tar.gz',
  'pd.rta.1.0_3.12.2.tar.gz',
  'pd.rusgene.1.0.st_3.12.0.tar.gz',
  'pd.rusgene.1.1.st_3.12.0.tar.gz',
  'pd.s.aureus_3.12.0.tar.gz',
  'pd.soybean_3.12.0.tar.gz',
  'pd.soygene.1.0.st_3.12.0.tar.gz',
  'pd.soygene.1.1.st_3.12.0.tar.gz',
  'pd.sugar.cane_3.12.0.tar.gz',
  'pd.tomato_3.12.0.tar.gz',
  'pd.u133.x3p_3.12.0.tar.gz',
  'pd.vitis.vinifera_3.12.0.tar.gz',
  'pd.wheat_3.12.0.tar.gz',
  'pd.x.laevis.2_3.12.0.tar.gz',
  'pd.x.tropicalis_3.12.0.tar.gz',
  'pd.xenopus.laevis_3.12.0.tar.gz',
  'pd.yeast.2_3.12.0.tar.gz',
  'pd.yg.s98_3.12.0.tar.gz',
  'pd.zebgene.1.0.st_3.12.0.tar.gz',
  'pd.zebgene.1.1.st_3.12.0.tar.gz',
  'pd.zebrafish_3.12.0.tar.gz'
)
install_with_url(annotation_url, pd_annotation_pkgs)

# Load this libraries because apparently just installing it isn't
# enough to verify that the correct versions of dependencies are installed.
library('foreach')
