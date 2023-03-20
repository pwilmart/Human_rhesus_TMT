# Human_rhesus_TMT

Analysis discussion of a multi-sample, multi-fraction, multi-kit, multi-species TMTpro experiment. How to analyze a 21 rhesus samples, 24 human samples, 45 samples total, 17 channels per plex (15 plus 2 pooled standards) in 3 plexes labeled with TMTpro 18-plex reagents experiment.

## Introduction

Animal models are common in human health research for many reasons that we will not go into. There are ethical arguments for and against. We will not explore that, either.

This core project involved amniotic fluid (AF) samples from rhesus monkeys and from some matched human subjects. There were three gestational ages when AF was collected. There were 21 rhesus samples and 24 human samples. The questions were how does gestational age alter the AF proteome for each species and how do the rhesus and human sample compare at each gestational age.

> These monkeys are part of a larger study to see effects of drug use on reproductive health. They have the sex and drugs parts well in hand. If they can form a 70s rock cover band, they will be living the dream...

We mostly do comparative proteomic experiments using TMT labeling and the SPS-MS3 acquisition on our Orbitrap Tribrids. We have more than 18 samples for each species, so we need to use multiple TMT kits (referred to as "plexes"). We could do the rhesus in two plexes and the human in two plexes. That would be 4 plexes total. We would have to figure out some way to compare the rhesus results to the human results from two independent multi-plex TMT experiments.

In the past, we combined rhesus and human samples into single plex TMT experiments (when we had fewer samples of each). Those could be analyzed by processing the data twice. Once using a human FASTA file (for the human sample channels) and a second time with a rhesus FASTA file (for the rhesus sample channels). This reduced TMT reagent cost and instrument time compared to doing each species in separate TMT experiments. Could we do something similar for this 45 sample experiment?

We could fit all samples into 3 TMTpro 18-plex kits (15 samples per plex with 2 pooled standard channels). This would reduce costs of reagents and instrument time compared to 4 plexes (independent processing of samples from each species would require 2 plexes each). The only real difference to previous single plex experiments would be having a pooled standard for internal reference scaling (IRS) made up of both human and rhesus protein digests. The pooled standard would have all possible AF peptides from both species, so this seems possible.

> Another relevant thought about combining rhesus and human samples in the same plexes comes from the nature of TMT experiments. TMT has high precision for data measured in the same scan. The data collected for the peptides that are identical between rhesus and human will be data that comes from the same scans. It is potentially more precise data. Note that the intensities from the species-identical peptides will be a subset of the total data, though.

## Data collection

Amniotic fluid is another wide dynamic range bio fluid with a lot of serum albumin. Serum albumin was depleted in all 45 samples using a human albumin depletion kit (we verified that it worked fine for rhesus, too). Determining how much of each sample to mix together in each plex was more complicated. We typically do a single shot (one LC run) analysis of an equal volume mixing then determine the total reporter ion signal for all peptides per sample. We want equal total peptide digest, so we adjust mixing volumes based on the total intensities (assuming summed reporter ions is a good proxy for total peptide).

We get different total intensities by species depending on which species FASTA file we use in the analysis. Clearly, we want to match the human samples using the intensities from a human FASTA search. We also want to match the rhesus samples using the rhesus FASTA results. Matching the human and rhesus is where it got more arbitrary. The average rhesus intensity was a little lower (10-20%) than the average human intensity, so we also boosted all the rhesus samples to (hopefully) make them more similar to human levels.  

After figuring out the mixing for the 17 samples (15 biological samples and 2 pooled standards) per plex, about 50 micrograms of total digest per plex was fractionated into 20 high pH reverse fractions (an online 2D LC setup). There were about 350K MS2/MS3 scans acquired per plex on an Orbitrap Fusion Tribrid. Even with the serum albumin depletion, the overall ID rates were low (about 14% for human searches and 11% for rhesus searches). These ID rates were similar to what we have seen in previous AF experiments and to what we see in other bio fluids like plasma. The canonical rhesus reference "one protein per gene" FASTA file has a similar number of protein sequences as human, but we consistently see significantly fewer PSMs identified at 1% FDR for rhesus searches compared to human searches (remembering to keep the human and rhesus samples in mind when counting things). My suspicion is that the quality of the human FASTA file (or choices of canonical sequences) is better for human (at least for the proteins present in AF).

## Analysis plan

Comparisons of the AF proteomes across gestational ages for each species is straightforward. We get the total protein reporter ion intensities for each species in separate analyses using the respective species FASTA file. We do the same basic IRS workflow for each FASTA search results:

- pick a minimally redundant (from tryptic peptide point of view) FASTA file that has high genome coverage. These are typically the "one protein per gene" reference proteomes with canonical sequences (no protein isoforms).

- Do Comet searches with wide parent ion tolerance and ion trap fragment ion settings. Use static mods for TMT reagents and only variable oxidized Met. Use semi-tryptic cleavage because this is a bio fluid.

- Filter PSMs at 1% FDR using delta mass conditional score distributions and the target/decoy method (concatenated, reversed protein sequences).

- Infer proteins using basic parsimony protein grouping and extended parsimony protein grouping. Define unique (usable) and shared peptides (not usable) for quantification based on the final inferred list of proteins, protein groups, or protein families.

- Sum reporter ion intensities from all scans associated with usable peptides per protein per channel.

- Perform normalizations, IRS methods, and statistical testing on protein-level results.  

> We do not want use a two-species concatenated FASTA file (at least not with my pipeline). That would have all the identical peptides between species mapping to both a human and a rhesus protein entry. In most cases, those proteins would have peptides unique to each species as well. Any distinguishable proteins would have many shared peptides (shared between species) that would be excluded from quantitative analysis. We would lose a large fraction of the quantitative data.   

When we do the rhesus FASTA search, we get larger values for the rhesus channel intensity totals. Along the same lines, when we use the human FASTA file in the processing, we get larger values for the human channels. About 2/3 of the total reporter ion signal comes from peptides that are identical between species. The remaining 1/3 comes from the peptides that are different in each species. That is why we see larger total intensities for rhesus samples with a rhesus FASTA file (where we pick up the unique-to-rhesus peptides) and larger total intensities for human samples with a human FASTA.

Even though the pooled standards have different intensity values depending on the FASTA file used, the standards are the same in all plexes and the values of the reference channels still put all plexes on the same grand intensity scale **for each respective FASTA file search**. We search the 3 plexes once with a rhesus FASTA file, do the IRS adjustment using the pooled standard channels, and then compare the 21 rhesus samples by gestational age. Human samples are ignored. We do a second analysis using a human FASTA file, do the IRS adjustment using different pooled standard channel intensities, and then compare the 24 human samples by gestational age. Rhesus samples are ignored in this second analysis.

We work up the data from each FASTA search results using the same QC notebook and edgeR statistical testing notebooks used in most core projects. The notebooks were qualitatively very similar for either species and looked like the typical notebooks for projects where we had only one species. It is possible to multiplex species in large scale TMT experiments. While it works, it is probably getting too complicated. It is probably better to keep species separate, if possible.

> Some experiments have multiple species that you cannot process separately. Environmental metaproteomic samples are one example. Cancer xenograft systems are another example. Both of these examples are more complicated to analyze than a mixture of rhesus and human samples.  

## Separate Human and Rhesus Analyses

Numbers of PSMs directly influence the number of inferred proteins. The PAW pipeline uses the two peptide rule because that actually works. We had about 1450 proteins per plex for human searches. There were over 80 immunoglobulin proteins matched (about 5% of the total IDs). This is an artifact of how immunoglobulins are represented in the FASTA files. This needs to be kept in mind anytime you are looking at secreted bio fluids. They all have immunoglobulins present and the immunoglobulin protein matches will inflate protein ID numbers and complicate quantitative comparisons. After running the IRS procedure on the 3 plexes, there were just under 1400 quantifiable human proteins. The typical QC metrics were applied to the 24 human samples from the human FASTA searches and everything looked like a typical single-species IRS experiment. We can safely extract the human sample data from the two-specie IRS experiment.

The rhesus searches had fewer filtered PSMs compared to human and that translated into fewer protein IDs. The average number of IDs per plex was a little over 1300. There were just over 60 immunoglobulins for the rhesus samples. After IRS, there were just under 1300 quantifiable proteins. Again, the QC metrics for the 21 rhesus samples were like any IRS experiment. The statistical testing notebook indicated that the final IRS-adjusted rhesus data resembled other single-species experiments. Combining two species in the same TMT experiment did not seem to have any negative consequences when two independent analyses (one for each species) of the data were done.   

> Proteomics spent most of its formative years preoccupied with identifying things and most data analysis software reflects that. An ID first mentality is not a good approach for quantitative proteomics. In the human samples, we had over 1800 total identifications from the 3 plexes (a union of all protein IDs using the conservative 2 peptide rule). We averaged about 1450 protein IDs per plex. This average sample ID number is a good upper bound to what might be the number of quantifiable proteins in the experiment. We will always have more things we can ID than we can quantify (quantification is a higher bar). After the IRS procedure on the 3-plex human results, there were just under 1400 quantifiable proteins. <br><br> The emphasis on the union of protein IDs as the main metric of an experiment leads to thinking that quantification of 1400 proteins out of 1800 is not very good (and maybe that shiny new method is the answer), rather than thinking 1400 out of 1450 was **amazing**. Instead of managing quantitative expectations, we manage to create disappointment...

## Human to rhesus comparisons

There are at least two options to explore for comparison of the rhesus samples to the human samples at the three gestational ages. We have many identical tryptic peptides between species (about 60-70% of the total intensity). We could use the summed reporter ions from the cross-species shared peptides (and the cross-species protein matching that they provide) to get protein intensity subtotals that may be directly comparable. That would require a custom version of the TMT intensity rollup script but seems possible to try (this is ongoing).

A more global proteome comparison can also be done from the two separate FASTA search results where we have independent protein IDs and quantitative values. We need some way to align the rhesus protein IDs to the human protein IDs so that we can compare the protein intensities. We can BLAST rhesus protein sequences against human protein sequences and try to find ortholog pairs (reciprocal best matches). There are scripts at https://github.com/pwilmart/PAW_BLAST that will BLAST one list of proteins against another list of proteins and find matches. This is done by making subset FASTA sequence files from each protein list and calling BLASTP with the subset FASTA files as arguments.

The next question is what protein lists do we use for each species? We could compare all rhesus canonical gene products (21K sequences) to all human canonical gene products (21K sequences). We would need to filter the matches down to the rhesus AF proteins that we identified, find the human match, then see if the human protein was in the human AF list.

This does not fully leverage all the work we did to get the lists of identified amniotic fluid proteins for each species. Amniotic fluid has unique functional properties related to the specific proteins in AF. It seems logical that this special function would be maintained between species by having similar proteins in the proteome. If we restrict the two protein lists to the identified AF proteins, we simplify the task of deciding which AF proteins match and which do not match.

We really have two lists of AF proteins for each species. We have all the proteins (and protein groups) that we could confidently identify (the union of IDs across the three plexes). Those lists are 1800+ and about 1650 for human and rhesus, respectively. We also have the smaller lists of protein that we could quantify (just under 1400 and 1300, respectively). To make a long story shorter, it works better to define the ortholog relationships using the longer lists of identifications and then filter down to the quantifiable proteins.

We made subset FASTA files with AF IDs for human and rhesus, respectively. We picked rhesus as the reference proteome for a couple of reasons. We want to know how similar rhesus AF is to human AF for use as an animal model for human pregnancy health studies. Comparisons are also cleaner when comparing a shorter list to a longer list. We did the BLAST runs with the rhesus sequences as the query and the human sequences as the hit choices (query and hit are terms that BLAST uses). This creates an ortholog map from rhesus to human AF proteins (nearly all rhesus AF proteins have strong BLAST alignments to human AF counterparts). We have another script (https://github.com/pwilmart/annotations) that can add the ortholog information to the rhesus results and then add extra annotations from the human orthologs to the rhesus results table.

Starting with a copy of rhesus results, we can delete columns that we do not need to simplify the results table to some rhesus protein ID information (accessions and descriptions), the TMM-normalized IRS-adjusted reporter ions for each of the 21 samples, and human ortholog mapping information. We create a similar simplified results table from the human FASTA search results and add that as a tab to the sheet. We will use this data and the VLOOKUP function to align rhesus protein rows to their human ortholog rows.

We need to do some more data wrangling before we merge the rhesus and human tables. We have many immunoglobulin protein matches in both species (83 for human and 65 for rhesus) that we do not want to align individually (it will be error prone). We will get more reliable merging if we group all immunoglobulins together (by summing the intensities) and have one immunoglobulin group per species.

The ortholog mapping script is designed to facilitate adding protein annotation information to species with poorly annotated sequences. Multiple rhesus proteins can map to the same human ortholog. For the handful of human proteins that were present in the ortholog map list multiple times, the matches were reconciled to a single one-to-one best match. This moved a few rhesus proteins from having a human ortholog match to not having a match.

Now the list of about 1200 quantifiable rhesus proteins (with mmunoglobulins replace by a single group) could be used to see if the human ortholog had a match to the list of about 1300 quantifiable human AF proteins (with immunoglobulins replaced by a single group) using the VLOOKUP function in Excel. The lookup value to add to the rhesus table was a concatenated text string of the human reporter ions for all 24 human samples for each protein. After adding the matching values to the rhesus table, the human data was expanded back to 24 columns of reporter ion quantities.

> We now have data rows with 21 rhesus protein total intensities and the 24 human protein total intensities from the ortholog pair (nearly all rhesus protein had human orthologs). These values are not directly comparable between species for something like differential abundance testing. The rhesus proteins are sums of two kinds of peptide signals - those shared with human proteins (about 60% on average) plus all unique to rhesus peptides for that protein. Similarly, the human total intensities are sums of the cross-species shared peptide intensities and unique-to-human peptide intensities. Because of the respective unique peptide totals, we may end up with rather different protein intensity totals. This cannot be remedied by conventional normalization methods. It is something like the plex effect that we could handle with pooled reference channels. We cannot make a magic sample like the pooled standard that would be “the same” between two different species. If we happened to have identical proteomes for rhesus AF and human AF (all proteins being ortholog pairs with the same relative abundances), we would measure something that would look like IRS reference channels *before* the IRS adjustment (something with a lot of scatter in a scatter plot). We also expect that the two proteomes would not be perfectly similar and that would add more scatter.

If we cannot answer questions about which rhesus and human protein are at similar or different relative abundances, what do we do? What are the questions we can answer with this ortholog-matched quantitative data? We need to focus on global questions. We can ask questions like what fraction of the rhesus AF proteome had quantifiable human AF orthologs. And the related question of what fraction of the human AF proteome had rhesus AF orthologs. We can explore those questions for each of the three gestational ages. We can average intensities over all samples at each gestational age for each species and make scatter plots (by gestational age) of the average intensities. This might reduce some of the variability and let us see if the proteome similarity changes with gestational age.

The questions of what fraction of the respective proteomes overlap with the other species depends on how you count things. The most common approach is to count protein IDs (each protein weighted by 1) in some sort of Venn diagram presentation. We have a quantitative dimension to exploit. We can weight each protein by its average total intensity and sum those by overlap out of the total intensity. These intensity values are different for each species, so we do not get the same fractions in the overlap. We need to formulate questions like what fraction of the rhesus AF proteome overlaps with the human AF proteome (by intensity), and what fraction of the human AF proteome overlaps with the rhesus AF proteome (by intensity).

Another reason to think about using intensities when we are doing TMT is that the biological samples (and the gestational ages) are hidden from the peptide ID and protein inference. We can only compare protein ID overlap at the experiment level (where individual samples are not distinguishable). When we use reporter ion intensity totals for overlap, we can compute things for each sample and do things like averaging over gestational ages.

Protein ID counts are like spectral counting in that they underestimate high abundance and overestimate low abundance. For example, 95% of the rhesus proteins have human ortholog pairs. In the other direction, 88% of the human proteins have rhesus ortholog pairs. This also illustrated how list length affects apparent overlap. When we use total reporter ion totals, the intensity overlap is around 99% (and much more symmetric between comparison direction). There is a high degree of overlap in the AF protein lists and abundances between rhesus and human.

I am still working on trying to directly compare the abundance levels of ortholog pairs between the species. I do not yet have any prediction if that will be possible or not. At this point in time, I have not ruled it out.

## Data Analysis Steps (so far)

A few months ago I did a Twitter poll about factors that prevented folks from doing more TMT labeling. The high cost of the experiment was the clear winner in the poll. At the bottom was how difficult it is to analyze TMT experiments. Since I spend all my time analyzing TMT experiments, and I know how much of my time and energy that takes, I was shocked. Maybe my idea of what analyzing data means is different from the community norm. That would not be surprising because almost nothing I do follows community norms. Here is a very much abbreviated list of the data analysis steps on this project (and I am not done yet):

1. Download RAW files for each plex (3 folders)
2. Convert RAW files to MSn files for each plex (3 folders)
3. Duplicate MSn files for each plex (6 folders)
4. Denote one set of MSn files for human searches (3 folders)
5. Denote one set of MSn files for rhesus searches (3 folders)
6. Perform human DB semi-tryptic Comet searches for human (3 folders)
7. Perform rhesus DB semi-tryptic Comet searches for rhesus (3 folders)
8. Filter human DB semi-tryptic Comet search results for human (3 folders)
9. Filter rhesus DB semi-tryptic Comet searches for rhesus (3 folders)
10. Infer proteins for human DB semi-tryptic Comet search results for human (3 folders)
11. Infer proteins for rhesus DB semi-tryptic Comet searches for rhesus (3 folders)
12. Compile dataset stats for human DB track
13. Compile dataset stats for rhesus DB track
14. Copy all filtered human files (3 folders) to combined human folder
15. Copy all filtered rhesus files (3 folders) to combined rhesus folder
16. Infer proteins from combined human searches (1 folder)
17. Add TMT intensities to human protein summary
18. Infer proteins from combined rhesus searches (1 folder)
19. Add TMT intensities to rhesus protein summary
20. Make copy of human protein summary table and add sample labels to TMT channels
21. Flag any additional contaminant proteins in human summary
22. Zero out any IRS channels for human proteins with a missing IRS value
23. Run IRS script on human results
24. Add additional missing intensity protein labels to human IRS results files
25. Make copy of rhesus protein summary table and add sample labels to TMT channels
26. Flag any additional contaminant proteins in rhesus summary
27. Zero out any IRS channels for rhesus proteins with a missing IRS value
28. Run IRS script on rhesus results
29. Add additional missing intensity protein labels to rhesus IRS results files
30. Create and edit human QC notebook for IRS results
31. Create and edit rhesus QC notebook for IRS results
32. Create and edit human edgeR exact test for gestational times
33. Create and edit rhesus edgeR exact test for gestational times
34. Compile DE candidate stats by gestational time and species
35. Add DE candidate stats to dataset stats workbook
36. Create and edit human edgeR glm ANOVA test for gestational times
37. Create and edit rhesus edgeR glm ANOVA test for gestational times
38. Create Excel file of IRS output for human FASTA results and get intensity totals by species
39. Create Excel file of IRS output for rhesus FASTA results and get intensity totals by species
40. Add PSM/peptide intensity analysis tab to data set stats workbook
41. Get list of human protein IDs from human FASTA searches and make human FASTA subset
42. Get list of rhesus protein IDs from rehesus FASTA searches and make rhesus FASTA subset
43. BLAST human subset FASTA versus rhesus subset FASTA to find orthologs
44. BLAST rhesus subset FASTA versus human subset FASTA to find orthologs
45. Select some proteomics columns from human results for left columns in edgeR results sheet
46. Add edgeR exact test results to sheet and check row alignment
47. Add edgeR glm ANOVA test results to sheet and check row alignment
48. Format sheet to increase readability and add some columns for the ANOVA results
49. Add extra annotation columns to sheet
50. Add ReadMe column key tab
51. Select some proteomics columns from rhesus results for left columns in edgeR results sheet
52. Add edgeR exact test results to sheet and check row alignment
53. Add edgeR glm ANOVA test results to sheet and check row alignment
54. Format sheet to increase readability and add some columns for the ANOVA results
55. Add extra annotation columns to sheet for human orthologs
56. Add ReadMe column key tab
57. Start with copy of rhesus results
58. Simplify table to basic protein info, human ortholog, and IRS quant data
59. Add tab with human IRS quant data and human accessions
60. Find/flag human immunoglobulins (including IgA) and combine intensities
61. Replace individual human Ig proteins with combined total
62. Find and replace rhesus Ig proteins with combined rhesus total
63. Make key, value columns for the human ortholog data in the human tab
64. Add VLOOKUP column to rhesus tab to fetch data for human orthologs
65. Expand the lookup values back into human quant columns
66. Find redundant ortholog mappings and reconcile to one-to-one mappings
67. Compile some intensity based statistics on what fraction of proteomes are mapped
68. Make a scatter plot of protein intensities for matched ortholog pairs for G85 age
69. Average intensities for rhesus and for human for all G85 samples, respectively
70. Copy data to scatter plot template
71. Add human gene symbol column
72. Need to exclude top protein in rhesus (very different abundance in human)
73. Need to match median intensities so that slope is closer to 1
74. Sort table by decreasing rhesus intensity
75. Label axes and add plot title
76. Make similar scatter plot for G110 ages
77. Make similar scatter plot for G135 ages
