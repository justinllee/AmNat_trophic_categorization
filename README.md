# Trophic niche representation affects macroevolutionary inference 

The follwing is a guide to the supplemental material for the manuscript entitled "Trophic niche representation affects macroevolutionary inference: a case study in elapid snakes", which was accepted in the journal American Naturalist.



## Files included in Datasets.zip :

- **'All_Elapids_squamatabase_JLL.csv':**
Dataset containing all 7,890 predator-prey observations of snakes used in study, including dissected stomach contents (noted by "this study" in references tab), and primary literature data. All predator-prey observations are notated in SquamataBase format (see Grundler 2020). Seven additional columns are included that note the diet category schemes used throughout the study, and the subgroup (noted as 'radiation') of each elapid. Explanations of all database fields/columns in the dataset are summarized below (refer to the SquamataBase GitHub page: https://github.com/blueraleigh/squamatabase, for more details):

- **'All_Elapids_tree_v2.tre':**
Modified phylogenetic tree of Elapidae based on the phylogeny of Title et al. (2024) (original tree is included as file 

- **'best_ultrametric_fulltree_ddBD_revision.tre':**
This tree also includes enhancements made based on robust evidence from recent taxonomic and systematic literature. Details of enhancements are included in the Materials and Methods section of the manuscript text, and code used to enhance the tree are included in the R-script 

- **'Script1_correct_Elapid_tree.R':**
A full list of species included in the phylogeny can be found in file  'Elapid_taxon_list_intree.csv'. 

- **'best_ultrametric_fulltree_ddBD_revision.tre':**
Raw phylogenetic tree of squamates from Title et al. (2024) (see references section in manuscript text).

- **'Elapid_full_taxon_list.csv':**
Entire list of elapid snake species included in phylogeny ('All_Elapids_tree_v2.tre') or in the dietary dataset ('All_Elapids_squamatabase.csv'). Original phylogeny from Title et al. (2024) (see references section in manuscript text).

- **'Elapid_taxon_list_intree.csv':**
List of all 245 elapid snake species included in the enhanced phylogeny (see file 'All_Elapids_tree_v2.tre') published by Title et al. (2024) (see references section in manuscript text).



## R-scripts:

- **'Script1_Correct_elapid_tree.R':**
This script enhances the Title et al. (2024) phylogeny so that recent taxonomic and systematic changes are included in the Elapid tree.

- **'Script2_Generate_diet_datasets.R':**
This script merges the Elapid predator-prey observations and diet categories from the file "All_Elapids_squamatabase.csv" into six separate datasets. They are used to calculate diet niche states using the Dirichlet Process Multinomial State (DPMS) model (Grundler and Rabosky 2020). See scripts 'Script3_Figure2_Phylo_DPMS_analysis.R' and 'Script12_FigureS4_Nonphylo_DPMS_analysis.R'.

- **'Script3_Phylo_DPMS_analysis.R':**
Script to estimate dietary niche states of elapids using phylogenetic implementation of DPMS model.

- **'Script4_Figure2_Calculate_diet_breadth.R':**
Script to estimate dietary breadth scores based on phylogenetic implementation of DPMS model used in Figure 2.

- **'Script5_Figure3_Generate_rarefaction_curves.R':**
Script to generate rarefaction curves to compare diets of hydrophiines and coralsnakes used in Figure 3.

- **'Script6_Figure4_Calculate_diet_tiprates.R':**
Script to estimate diet tip rates of Elapids based on DPMS analyses in Figure 4 and diet breadth scores (plots for diet breadth are coded separately; see Figure S6). Plot tip rates across phylogeny of Elapidae. Comparisons between diet breadth/proportions and their rates are also calculated.

- **'Script7_FigureS1_(A)_Generate_category_T2_order.R':**
Generate phylogenetic [T2] category scheme based on order-rank data from Timetree of life and existing elapid diet data (Kumar et al. 2022). 

- **'Script8_FigureS1_(B)_Generate_category_T3_family.R':**
Generate phylogenetic [T3] category scheme based on family-rank data from Timetree of life (Kumar et al. 2022) and existing elapid diet data.

- **'Script9_FigureS2_(A)_Barplots_of_diet_states_phylo.R':**
Visualize dietary niche states estimated from phylogenetic implementation of DPMS model (Figure S2, in part).

- **'Script10_FigureS2_(B)_Barplots_of_diet_states_nonphylo.R':**
Visualize dietary niche states estimated from non-phylogenetic implementation of DPMS model (Figure S2, in part).

- **'Script11_FigureS3_Plot_Phylo_DPMS_analysis':**
Plot phylogenetic implementation of DPMS analysis onto a phylogeny of Elapidae.

- **'Script12_FigureS4_Nonphylo_DPMS_analysis.R':**
Estimate dietary niche states based on the non-phylogenetic implementation of DPMS model. Script also generates non-phylogenetic dietary breadth scores based on the model and plots dietary niche states onto a phylogeny of Elapidae.

- **'Script13_FigureS5_Nonphylo_diet_breadth.R':**
Generate boxplots of dietary breadth scores estimated using non-phylogenetic implementation of Dirichlet-multinomial model. 

- **'Script14_FigureS6_Plot_diet_breadth_tip_rates.R':**
Plot tip rates calculated from dietary breadth scores in previous script ('Script13_FigureS5_Nonphylo_diet_breadth.R').



## R-packages (dependencies not listed):
- phytools_2.4-4
- ape_5.8-1 
- macroevolution_1.0
- RColorBrewer_1.1-3
- devtools_2.4.5
- remotes_2.5.0
- scales_1.4.0
- classInt_0.4-11
- bm_1.0
- phylo_1.0

R version 4.5.0 (2025-04-11)

Running under: macOS Sequoia 15.7.7



## Supplemental Figures: 
**Figure S1:** Phylogenetic reconstruction of categories at order rank [T2] (A) and family 
rank [T3] (B). Clades younger than the K-Pg extinction were grouped together and given the 
closest clade-level name corresponding to them.

**Figure S2:** Barplots depicting elapid prey category proportions of inferred dietary niche 
states. Niche states were inferred under a phylogenetic (A) and non-phylogenetic (B) 
DPMS model for six prey categorization schemes. All color legends correspond to prey 
categories except categories T2–T3, which were too large to depict  (gradients were used 
in their place and correspond to the same colors used in category T1). 

**Figure S3:** Phylogeny of elapids with dietary states inferred under six category schemes. 
Each ring of shapes correspond to niche states inferred by the phylogenetic implementation 
of the Dirichlet process multinomial state model (DPMS) model under different 
categorizations (categories from inner-to-outer ring: T1, T2, T3, E1, E2, E3; see bottom 
right legend). Colors at each ring of tips indicate different dietary states inferred for 
each category. Colors on tree branches indicate the six elapid subgroups grouped in this 
study (Pink = Afro-Asian elapids; Yellow = Asian coralsnakes; Orange = Australo-Papuan 
hydrophiines; Purple = kraits; Turquoise = New World coralsnakes; Green = sea snakes; see 
bottom left legend). Refer to figure S4 for results under the non-phylogenetic DPMS 
analysis.

**Figure S4:** Phylogeny of Elapidae with dietary niche states inferred by the 
non-phylogenetic implementation of the DPMS model. Colors and shapes correspond with 
Figure S3.

**Figure S5:** Box and whisker plots showing differences in log transformed dietary breadth 
amongst elapid radiations based on a non-phylogenetic Dirichlet-multinomial model. Points 
represent elapid species. Refer to Figure 3 for colors and letters.

**Figure S6:** Rates of diet evolution amongst Elapidae based on diet tip rates of all six 
prey categorization schemes. Data used to generate tip rates were the raw diet breadth for 
each species.



## References:
Grundler, M.C. 2020. SquamataBase: a natural history database and R package for 
comparative biology of snake feeding habits. Biodiversity Data Journal 8:e49943.
https://doi.org/10.3897/BDJ.8.e49943. 
[GitHub repository: https://github.com/blueraleigh/squamatabase]

Grundler, M.C., and D.L. Rabosky. 2020. Complex ecological phenotypes on phylogenetic 
trees: a Markov process model for comparative analysis of multivariate count data. 
Systematic Biology 69:1200–1211. https://doi.org/10.1093/sysbio/syaa031.
[GitHub repository: https://github.com/blueraleigh/macroevolution]

Kumar, S., M. Suleski, J.M. Craig, A.E. Kasprowicz, M. Sanderford, M. Li, G. Stecher, and 
S.B. Hedges. 2022. TimeTree 5: an expanded resource for species divergence times. 
Molecular Biology and Evolution 39:msac174. https://doi.org/10.1093/molbev/msac174. 
[https://timetree.org/]

Title, P.O., S. Singhal, M.C. Grundler, G.C. Costa, R.A. Pyron, T.J. Colston,, M.R. 
Grundler, et al. 2024. The macroevolutionary singularity of snakes. Science 383:918–923.
[GitHub repository: https://github.com/macroevolution/squamata]



## Appendix:

#### Database fields used in 'All_Elapids_squamatabase_JLL.csv'
  - no:

    The number of existing diet records in SquamataBase. New diet records are noted as NA.
  
  - predator:
  
    Scientific name of predator. Some names have been updated to match the tip labels of the Title et al. (2024) tree, or to match recent taxonomic updates.
  
  - predator_verbatim:
  
    Scientific name of predator reported in the original source of the diet record. 
  
  - predator_rank:
  
    Linnean rank of predator.
  
  - predator_taxon:
  
    Semicolon separated list of taxonomic ranks applicable to predator. 	
  
  - predator_age:
  
    Age of predator specimen.
  
  - predator_sex:
  
    Sex of predator specimen.
  
  - predator_count:
  
    Number of individual predators involved in diet record.
  
  - predator_mass:
  
    Mass (in grams) of predator.
  
  - predator_svl:
  
    Snout-to-vent length (generally, from snout tip to vent) of predator. All measurements in millimeters.
  
  - predator_tl:
  
    Total length (snout-to-vent length + tail length) of predator. All  measurements in millimeters.
  
  - predator_voucher:
  
    A museum voucher number or random alphanumeric code applied for each predator specimen. Random alphanumeric codes were not given for new literature diet records reported herein.
  
  - prey:
  
    Scientific name of predator. Some names have been updated to match recent taxonomic updates.
  
  - prey_verbatim:
  
    Scientific name of prey reported in the original source of the diet record. 
  
  - prey_rank:

    Linnean rank of prey.
  
  - prey_taxon:

    Semicolon separated list of taxonomic ranks applicable to prey.
  
  - prey_age:

    Age of individual prey item.
  
  - prey_sex:

    Sex of individual prey item.
  
  - prey_count:
  
    Number of prey items involved in diet record.
  
  - prey_mass:

    Mass (in grams) of prey item.
  
  - prey_tl:

    Total length of prey item. All measurements in millimeters.
  
  - prey_ingested:

    Direction prey item was ingested (e.g., head first, tail first, bent double).
  
  - prey_voucher:
  
    A museum voucher number or random alphanumeric code applied for each prey item. Random alphanumeric codes were not given for new literature diet records reported herein.
  
  - locality_adm0_name:
  
    Country where diet record occurred.
  
  - locality_adm1_name:

    State/provincial level region where diet record occurred.
  
  - locality_adm2_name:

    County/municipality level region where diet record occurred.
  
  - locality_misc:

    Additional information relevant to location of diet record.
  
  - locality_longitude:

    Longitude (in decimal form) of diet record.
  
  - locality_latitude:

    Latitude (in decimal form) of diet record.
  
  - event_basis:

    The basis for the diet record (e.g., direct observation, dissected stomach/gut contents, etc.)
  
  - event_setting:

    Note stating whether the diet record was observed naturally or in a captive setting.
  
  - event_date:

    The date when the diet record was presumed to have occurred (given in YYYY-MM-DD format). The collection date is given under this field for records derived from dissected stomach contents.
  
  - event_start:

    The time when the diet record started (in 24-hour HH:MM format).
  
  - event_end:

    The time when the diet record ended (in 24-hour HH:MM format).
  
  - event_outcome:

    The outcome of the predation event in the diet record (in most cases, this is prey_eaten, but direct observations based on field encounters may have other outcomes such as predation_interrupted_by_observer)
  
  - event_habitat:

    Habitat descriptor noting the setting of diet record (e.g., terrestrial, fossorial, arboreal, aquatic)
  
  - event_habitat_verbatim:

    Verbatim description of habitat for diet record.
  
  - event_remark:

    Any relevant notes of diet record that do not fall under other database fields.
  
  - reference:

    The reference/citation of diet record. Novel diet records reported by us are recorded using the entry this study.
  
  - prey_traditional:

    Dietary categorization of prey item under the traditional [T1] category scheme.
  
  - prey_phylo1:

    Dietary categorization of prey item under order rank [T2] scheme.
  
  - prey_phylo2:

    Dietary categorization of prey item under family rank [T3] scheme.
  
  - prey_functional1:

    Dietary categorization of prey item under functional/ecological [E1] category scheme.
  
  - prey_functional2:

    Dietary categorization of prey item under the habitat [E2] category scheme.
  
  - prey_type:

    Dietary categorization of prey item under the Mass-Bulk Theory (MBT) prey type [E3] category scheme.
  
  - radiation:

    Assigned elapid subgroup of each predator species for diet record.

#### Database fields and elapid shorthands used in 'Elapid_full_taxon_list.csv' and 'Elapid_taxon_list_intree.csv':

  - "taxon": Scientific name (at species rank) of elapid taxon.
	
  - "geog_radiation": Assigned elapid subgroup of each elapid species. 

  - 'afro-asian' = Afro-Asian Elapids (A)

  - 'asian-corals' = Asian Coral Snakes (B)

  - 'australo_papuan' = Australo-Papuan hydrophiines (C)

  - 'kraits' = Asian kraits (D)

  - 'nw-corals' = New World Coral Snakes (E)

  - 'sea_snakes' = Sea Snakes (F)
 
