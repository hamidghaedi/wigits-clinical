# Comparing our clonal haematopoiesis detection with the Vlasschaert approach

Prepared 26 November 2026 for the Rauh lab discussion. Worked on KQ0004, a glioma case with paired
blood and tumour whole-genome sequencing, blood mean depth 127.7x.

The reference method is Vlasschaert C, Mack T, Heimlich JB, Niroula A, Uddin MM, Weinstock J, Sharber B,
Silver AJ, Xu Y, Savona M, Gibson C, Lanktree MB, Rauh MJ, Ebert BL, Natarajan P, Jaiswal S, Bick AG.
*A practical approach to curate clonal hematopoiesis of indeterminate potential in human genetic
datasets.* Blood 2023;141(18):2214-2223, doi:10.1182/blood.2022018825, PMID 36652671. The supplementary
tables (BLOOD_BLD-2022-018825-mmc1.xlsx) were supplied by the Rauh lab on 26 November 2026 and are
incorporated in section 6 below, which supersedes the gene-level comparison in sections 1 to 5.

## 1. The two designs answer different questions

The Vlasschaert pipeline was built to ascertain clonal haematopoiesis across roughly 550,000 UK Biobank
and All of Us participants, where the output feeds epidemiological association testing. That purpose
demands one clean binary per participant, and a false positive is more costly than a false negative
because it dilutes effect sizes. The paper makes exactly this point, that small changes in filtering
parameters considerably increase misclassification and reduce the effect size of associations.

Our layer was built to report one patient at a time alongside a clinical genomic report. There the
reader is a person deciding what to do about an individual, a false negative is as costly as a false
positive, and the sensible output is a graded inventory with the evidence attached rather than a
binary. Neither design is a better version of the other.

| | Vlasschaert approach | This layer |
|---|---|---|
| Purpose | cohort ascertainment for association testing | per-case reporting |
| Output | a filtered call set | every call, with a confidence band |
| Caller | Mutect2, tumour-only on blood | SAGE 5.0.2, blood in the sample-of-interest slot |
| Annotation | ANNOVAR | PAVE with gnomAD and ClinVar |
| Target definition | curated list of CH driver mutations, variant level | 50-gene panel, gene level |
| Depth requirement | 20x or more | none; the floor is computed and reported |
| Supporting reads | 5 or more | none; 5 and 3 define bands |
| Allele fraction | 2% or more | none; 2% and 5% define bands |
| Strand requirement | at least one read in each orientation | SAGE's statistical strand-bias test |
| Repeat handling | homopolymer indels at 3 bp or longer removed | SAGE repeat context and REDUX jitter parameters annotated |
| Artefact control | recurrent artefact positions learned from 550,000 samples | per-variant quality flags only |
| Paired tumour | not available | used to annotate presence in the tumour |

## 2. What happens when their filters are applied to our calls

Their sequence-based criteria were applied in order to the 178 variants our caller produced.

| Step | Remaining | Removed |
|---|---|---|
| all calls | 178 | |
| depth 20x or more | 178 | 0 |
| alternative allele reads 5 or more | 77 | 101 |
| allele fraction 2% or more | 77 | 0 |
| homopolymer indels removed | 75 | 2 |
| restricted to clonal haematopoiesis genes | 48 | 27 |
| excluding variants consistent with an inherited origin | 9 | 39 |

Three things are worth drawing out.

**The five-read rule does nearly all the filtering and the 2% rule does none of it.** At this depth five
reads is already an allele fraction of 3.9%, so every variant that clears the read requirement clears
the allele fraction requirement automatically. The 2% threshold only becomes the binding constraint
above 250x. This is the same arithmetic as the power analysis: the nominal definition is not what the
filter enforces.

**Twenty-six clone candidates in panel genes are removed by the read requirement alone**, the largest
being NF1 p.Leu190* at 3.5%, SRCAP p.Asn2890His at 3.2%, CSF3R p.Val224Gly at 2.6%, CALR p.Phe74Cys at
2.6%, BCOR p.Lys1063Gln at 2.3% and DNMT3A p.Ile898Asn at 2.2%. Each sits above the 2% definition and
below five supporting reads. Under their pipeline these are silently absent; under ours they appear in
the supported band with a note that the evidence is below what this depth resolves confidently. For
cohort work the removal is the right call, because at four reads the false-positive rate is not
controllable across 550,000 samples. For a single patient the same variant is worth a reviewer's eye.

**The 27 variants removed by restricting to panel genes are our error, not a methodological difference.**
They sit in 16 genes, among them APC, BAP1, BARD1, BRCA1, BRCA2, BRIP1, MEN1, MLH1, PMS2, RB1, TSC2 and VHL, pulled in because the run used Hartwig's cancer-predisposition
hotspot list rather than a clonal haematopoiesis hotspot list.

## 3. The nine survivors, and why none of them should be believed

Applying their published criteria to our data yields nine clonal haematopoiesis calls.

| Gene | Change | Allele fraction | Alt reads | Depth | Also in tumour |
|---|---|---|---|---|---|
| RUNX1 | p.Ser424= | 0.187 | 15 | 82 | yes |
| CREBBP | p.Ala1944Gly | 0.131 | 14 | 106 | yes |
| CREBBP | p.Thr1932Pro | 0.098 | 12 | 122 | yes |
| CREBBP | p.Lys1495Asn | 0.095 | 14 | 150 | yes |
| SETBP1 | p.Leu1573Arg | 0.072 | 6 | 83 | yes |
| CEBPA | p.Phe33Cys | 0.060 | 7 | 117 | yes |
| NF1 | p.Met1? | 0.058 | 7 | 120 | yes |
| NF1 | p.Thr25Lys | 0.054 | 7 | 129 | yes |
| ATM | p.Ser421* | 0.042 | 6 | 143 | yes |

**All nine are present in the tumour.** A genuine blood clone is blood-enriched, so a variant matched
across both libraries at comparable fraction is more consistent with a systematic artefact than with
clonal haematopoiesis. Eight of the nine carry a low base-quality flag, four a poor mapping-quality
flag and two a strand-imbalance flag. One, NF1 p.Met1?, is a start-codon call, and one, RUNX1 p.Ser424=, is synonymous. Three CREBBP calls cluster in one
gene, which is the signature of a problematic region rather than three independent clones.

This is the clearest single result of the comparison. Their published criteria, applied faithfully,
return nine calls from this sample and the correct number is probably zero. The criteria are not at
fault; what is missing is the part of their method we could not apply, namely the curated variant list
and the recurrent-artefact positions learned from the large cohorts.

## 4. Where each approach is stronger

**Theirs.** Variant-level curation is the decisive advantage. A gene-level panel accepts any variant in
a panel gene, so a synonymous RUNX1 change and a DNMT3A R882 hotspot are treated alike; the nine
survivors above include a synonymous variant that no curated list would admit. Their recurrent-artefact
positions, learned across 550,000 samples, are something no single-case analysis can reproduce, and the
paper singles out TET2 and ASXL1 as needing exactly this. Their population-based filtering step has no
counterpart in our layer, which can only use per-variant gnomAD frequency.

**Ours.** The confidence banding keeps the 26 sub-threshold candidates visible with their evidence
rather than discarding them, which suits a per-patient report. The detection floor is computed from
each sample's own depth and printed on the report, so a negative result states what it means. The
strand test is a statistical one rather than a requirement of one read in each orientation, and it
flagged two variants that their orientation rule admits. SAGE brings duplicate marking, base-quality
recalibration and microsatellite jitter parameters from REDUX, which bear directly on indel calling in
repeats. And we have the paired tumour, which supplied the observation that sank all nine calls and
which connects directly to the infiltration question in Buttigieg et al., Cell Rep Med 2025;6(3):101989.

## 5. What to do about it

1. Obtain the curated gene and variant list from the Blood 2023 supplement and adopt it as the target
   definition. This is the largest single improvement available and removes most of section 3.
2. Replace the cancer-predisposition hotspot input with a clonal haematopoiesis hotspot list.
3. Implement the tumour-sharing rule. The column is populated on every row and no rule consumes it.
4. Obtain any lab-maintained recurrent-artefact positions for TET2 and ASXL1.
5. Run Mutect2 and ANNOVAR alongside SAGE on the same samples. The two callers have different error
   models, the comparison is cheap, and it settles how much of section 3 is caller-specific.

Item 5 is a self-contained first project for a new trainee and would produce a concordance table of the
kind this document only approximates from one side.


## 6. What the supplementary tables change

The supplement to the Blood 2023 paper was supplied on 26 November 2026. It contains six tables: the
curated list of candidate variants (S1), recurrent artefacts removed (S2), variants validated by their
association with age and the TERT risk allele rs7705526 (S3), variants failing their binomial test (S4),
a worked example of a low-quality sample (S5), and the genes they considered and rejected (S6). Having
them turns several of the judgements above from inference into measurement.

### 6.1 The target definition operates at the variant level, not the gene level

The most consequential correction. Section 1 recorded their target definition as a curated list of CH
driver mutations and ours as a 50-gene panel, but it understated the gap. Their list names 74 genes, and
for each gene it names the mutation classes that count in that gene. Forty-seven genes admit truncating
variants, forty-six carry an explicit list of missense positions, and many carry both. CEBPA and NF1
admit truncating variants only. SETBP1 admits seven missense positions and nothing else. CBL admits
missense only within the RING finger, residues 381 to 421. Synonymous variants are admitted nowhere.

A call therefore has to clear two independent requirements: its gene is on the list, and its specific
change is one the list recognises for that gene. Our layer applies only the first.

### 6.2 Panel overlap

Forty-two of our fifty genes are on their curated list. Eight are not: ATM, CALR, CHEK2, DDX41, SH2B3,
SRCAP, YLPM1 and ZBTB33. None of the eight appears in S6 either, so they were not tested and rejected,
they were outside the scope of the curation. Thirty-two of their genes are absent from our panel:
ASXL2, CBLB, CSF1R, CTCF, CUX1, EED, EP300, ETNK1, GATA1, GATA2, GATA3, GNA13, IKZF1, IKZF2, IKZF3, JAK1,
JAK3, KMT2A, KMT2D, PDS5B, PDSS2, PHIP, PRPF40B, PTEN, SETD2, SETDB1, SF1, SF3A1, SMC1A, STAG1, SUZ12 and
U2AF2. S6 also shows which candidate genes they rejected and why, with age-association odds ratios and
myeloid-cancer counts per cohort, so the exclusions are documented rather than assumed.

### 6.3 Applying the curated definition to KQ0004

Section 5 reported that their sequence filters reduced our 178 calls to nine clone candidates, all nine
present in the tumour, three of them clustered in CREBBP. Applying the variant-level definition as well:

| Route | Calls retained |
|---|---|
| Our bands, no filtering | 178 (117 clone candidates) |
| Their sequence filters alone | 9 |
| Their curated variant definition alone | 4 |
| Both | 0 |

Every one of the nine falls outside the curated definition, by five distinct routes. RUNX1 p.Ser424= is
synonymous. The three CREBBP missense calls and SETBP1 p.Leu1573Arg sit at positions the curated lists
for those genes do not name. CEBPA p.Phe33Cys and NF1 p.Thr25Lys are missense in genes that admit only
truncating variants. NF1 p.Met1? is a start-lost, which is not among the admitted classes anywhere.
ATM p.Ser421* is nonsense, but ATM is not a curated gene.

The four that do fall inside the definition are all below their read-count threshold: NF1 p.Leu190* at
3.5% with 4 supporting reads, TP53 p.Tyr220Ser at 1.5% with 2, NF1 p.Gly1051* at 0.8% with 1, and
PTEN p.Leu247* at 0.7% with 1. The last is an off-panel record that entered through the germline hotspot
file defect noted in section 4, and PTEN happens to be a curated gene.

The intersection is empty. The two routes are close to orthogonal on this case: the sequence filters keep
well-supported calls at positions the definition does not recognise, and the definition keeps recognised
changes that are too thinly supported to call. The honest reading of KQ0004 under their published method
is that there is no clonal haematopoiesis to report, which is what the earlier power analysis predicted
for a 128x blood genome.

*Corrected by section 7: NF1 p.Leu190* was also a hotspot-forced call and is not present after
the re-run.*

The nearest thing to a candidate is NF1 p.Leu190*, a nonsense change in a curated gene at 3.5% allele
fraction with 4 supporting reads, absent from the tumour. It misses their five-read rule by one read.
Our layer places it in the supported band. This single variant is the clearest illustration of where the
two designs part company, and it is worth putting to the Rauh lab directly.

### 6.4 The low-quality sample signature fires on KQ0004

*Corrected by section 7. The ATM observation below rests on calls that the wrong hotspot input had
forced, and does not survive the re-run.*

Table S5 shows a sample they classified as likely artefact: seven variants in SETD2, NF1, CUX1, EZH2 and
CREBBP, all between 4.9% and 6.2% allele fraction. The pattern they are describing is several
independent loss-of-function calls at closely similar allele fractions, which a library artefact
produces and several independent clones do not.

Implemented as a rule and run on KQ0004, the signature fires, and it fires on the panel records alone,
not only on the off-panel block. Our panel carries seven truncating calls, and five of them are ATM
nonsense variants between 0.9% and 4.2%: p.Glu299*, p.Ser421*, p.Gly1016*, p.Leu1457* and p.Trp2104*.
Five independent nonsense mutations in one gene in one blood sample is not a plausible clonal history.
ATM is one of the eight genes we added beyond their list, and it contributes five of our seven panel
truncating calls and one of the nine sequence-filter survivors. Their decision to leave ATM out is
supported by our own data.

### 6.5 What we are adopting

The curated list becomes the target definition, carried as a resource file rather than written into
code. Membership of S3 is recorded as positive evidence, since a variant that carries an age association
and a TERT association across two cohorts of a quarter of a million people is unlikely to be an
artefact. S2 and S4 become an artefact annotation. Consistent with the reporting rule in section 1 of
the proposal, none of these removes a call from the table: a call outside the curated definition is
reported with that stated, because a clinician reading one patient is entitled to see what was found
and why it was not counted.

The eight genes we added stay in the panel, each carrying a note that it is outside the curated
definition, with ATM additionally carrying the artefact observation above. The thirty-two genes we are
missing are added to the panel for the next run.

Three things are not yet done: the germline hotspot file is still the wrong input and still pulls in 103
off-panel records; the tumour-sharing column is populated but unused; and the whole comparison rests on
one glioma case at 128x, so the depth distribution of the bladder cohort has to be pulled before any of
these numbers are treated as general.


## 7. The re-run, and two corrections to section 6

The panel was widened to the 82 genes, the germline hotspot file was replaced, and SAGE and PAVE were
re-run on the same blood alignment. The replacement hotspot list keeps the positions that are both in
the CHIP panel and recognised by the curated definition: 1,475 sites across 29 genes, against the
36,614 cancer-predisposition sites that were being supplied before.

| | Before | After |
|---|---|---|
| Panel genes | 50 | 82 |
| Hotspot sites supplied to SAGE | 36,614 predisposition | 1,475 CHIP |
| Calls | 178 | 119 |
| Off-panel calls | 103 | 0 |
| Genes with a call | 28 of 50 | 52 of 82 |
| Clone candidates | 117 | 36 |
| Inside the curated definition | 4 | 3 |
| Low-quality sample signature | fires | does not fire |

### 7.1 The ATM pattern was our artefact, not ATM's

Section 6.4 reported five ATM nonsense calls between 0.9% and 4.2% and took them as evidence that the
curated list is right to leave ATM out. That reading was wrong. ATM is a cancer-predisposition gene, so
its positions were in the germline hotspot file, and SAGE evaluates every hotspot position at hotspot
sensitivity whatever the panel says. Four of the five calls were forced by that input. With the CHIP
hotspot list in place ATM has one call left, p.Ser421* at 4.2% with 6 supporting reads.

The sample-level artefact signature no longer fires. The whole sample now carries two truncating calls,
ATM p.Ser421* and KMT2D p.Leu804*, at 4.2% and 4.3%. We have no evidence either way on whether ATM
belongs in a CHIP panel, and the earlier claim is withdrawn.

### 7.2 NF1 p.Leu190* was also hotspot-forced

Section 6.3 named NF1 p.Leu190* as the nearest thing to a candidate on this case. NF1 is likewise a
predisposition gene, that position was likewise in the germline hotspot file, and the call does not
survive the re-run. The three NF1 calls that remain are p.Thr25Lys at 5.4% and two synonymous changes
at germline allele fractions.

This is worth stating plainly: a wrong resource file does not merely add noise at the edges, it
manufactures exactly the kind of call that looks most interesting, a truncating variant in a
recognisable gene at a plausible clonal allele fraction.

### 7.3 What the curated definition now admits

Three of the 119 calls are inside the definition, and all three are weak:

| Gene | Change | Blood | Tumour | SAGE |
|---|---|---|---|---|
| BRAF | p.Asp594Glu | 1.6%, 2 of 128 reads | 0.5% | minTumorQual, minTumorSupport |
| TP53 | p.Leu330Arg | 1.5%, 2 of 137 reads | 0.7% | minTumorQual, minTumorSupport |
| TP53 | p.Tyr220Ser | 1.5%, 2 of 135 reads | 1.9% | minMQF, minTumorQual, minTumorSupport |

Each rests on two supporting reads, each is in the provisional band, and each is present at a similar
low fraction in the tumour as well, which is what tumour material in the blood draw looks like. None
would pass their five-read rule. The conclusion of section 6.3 is unchanged and now rests on a cleaner
call set: KQ0004 has no clonal haematopoiesis to report.

### 7.4 One independent confirmation

CREBBP p.Ala1944Gly, one of the nine calls that passed their sequence filters, matches a recurrent
artefact in their table S2. Eight of those nine are still called after the re-run, all eight outside the
curated definition, and this one is documented as an artefact in two cohorts of a quarter of a million
people. That is the artefact blacklist doing the job it was built for on the first case we applied it to.
