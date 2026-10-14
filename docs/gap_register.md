# Gap register

Features of the Hartwig clinical workflow that are no longer public, were discontinued, or exist only in the private oncoact repository. Each entry records what is known, where that knowledge comes from, and what this project intends to do about it. Add rows as gaps are found; close rows by linking the commit or directory that fills them.

Status values. `open`, `in progress`, `done`, `wontfix`.

| ID | Gap | What is known | Source | Plan | Status |
|---|---|---|---|---|---|
| G01 | PROTECT evidence matching | Full matching rules, on-label logic, reporting caps and output schema documented in the 2022 README and source | hmftools tag `rose-v1.3`, `protect/` | Reimplement against the ORANGE JSON rather than raw tool outputs, keeping the rule set and TSV schema | open |
| G02 | SERVE knowledgebase harmonisation | Algorithm, supported sources, output files and transvar-based hotspot generation documented; the SERVE code at that tag is complete | hmftools tag `rose-v1.3`, `serve/` | Decide which sources are usable under our licences (CIViC and CGI are open; CKB, OncoKB and iClusion are not), then build a SERVE-compatible actionability directory from those | open |
| G03 | ROSE conclusion sentences | Schema of the rules table and matching logic documented; the table itself was never public | `rose/README.md` at `rose-v1.3` | Write our own rules table for the genes in the driver panel, starting from the README examples | open |
| G04 | Patient reporter PDF | Chapter order, inputs, consent and purity logic readable in the 2022 source; nothing after that | `patient-reporter/` at `rose-v1.3`, example report HMF-FOR-080 | Build a reporter that reads the ORANGE JSON plus PROTECT and ROSE output; use the published example report as the layout reference | open |
| G05 | Germline consent model | Three LIMS reporting levels and four per-gene conditions documented; gene table not public | `patient-reporter/germline/` | Define a consent level in the samplesheet and a versioned gene table in the repo | open |
| G06 | QC-fail report templates | Five fail reasons and their template codes documented | `patient-reporter/qcfail/`, OncoAct user manual | Implement fail reports as a separate path driven by a QC verdict file | open |
| G07 | Purity-unreliable handling | Rule is PURPLE fit method not NO_TUMOR; blanking behaviour documented | `patient-reporter/algo/QualityOverruleFunctions.java` | Reproduce the blanking rules; add the low-sensitivity template for 10 to 20 percent purity | open |
| G08 | Clinical trial matching | iClusion is a Dutch commercial feed; ACTIN is Hartwig's own trial study | `serve/README.md` | Out of scope for the Netherlands feeds; evaluate an open trial source later | open |
| G09 | SERVE known-event files feeding SAGE and LINX | SAGE hotspots and LINX known fusions in the public resource bundle are already derived from SERVE | `pipeline/README_RESOURCES.md`, OncoAct release notes 5.22 | Use the bundle as is; regenerate only if we change sources | open |
| G10 | ORANGE immune escape table | Computed into the JSON but the PDF rendering is commented out in ORANGE 5.0 | `orange/.../ImmunologyChapter.java` | Render it in our reporter | open |
| G11 | Cohort percentile files for ORANGE | Listed in the resource README but not read by current ORANGE code | `pipeline/README_RESOURCES.md`, `orange/` | Note only; no action unless a later ORANGE restores it | wontfix |
| G12 | Hartwig DOID resolution in ORANGE | DOID JSON loading is commented out, so configured DOIDs never resolve | `orange/.../OrangeAlgo.java` | PROTECT needs DOIDs, so resolve them in our layer from the samplesheet rather than relying on ORANGE | open |
| G13 | OncoAct version tracking | The clinical product has its own version line separate from WiGiTS | `pipeline5/ivdr/release-notes-5.22.md` | Keep a `VERSIONS.md` that pins oncoanalyser, every hmftools tool, and our layer | open |
| G14 | Report JSON written by the patient reporter | Existence known, schema unknown | `patient-reporter/` | Define our own schema; publish it with the reporter | open |
