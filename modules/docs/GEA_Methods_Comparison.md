# Comparison of Gene Expression Analysis Methods

This section compares the main RNA-seq quantification methods implemented in the pipeline, highlighting their workflows, strengths, and limitations.

---

## 1. HISAT2 Reference Guided

- **Workflow**: Aligns RNA-seq reads to a reference genome using HISAT2, then assembles and quantifies transcripts with StringTie.
- **Strengths**: High accuracy when a well-annotated reference genome is available; enables detection of known and novel transcripts.
- **Limitations**: Dependent on reference quality; may miss transcripts not present in the reference.

## 2. HISAT2 De Novo

- **Workflow**: Aligns reads to a de novo assembled transcriptome (no reference genome required), followed by quantification with StringTie.
- **Strengths**: Suitable for organisms lacking a reference genome; can discover novel transcripts.
- **Limitations**: Assembly quality depends on read depth and complexity; may produce fragmented transcripts.

## 3. Trinity De Novo

- **Workflow**: Uses Trinity to assemble transcripts de novo, then quantifies expression with StringTie.
- **Strengths**: Robust de novo assembly; widely used for transcript discovery in non-model organisms.
- **Limitations**: Computationally intensive; may generate redundant or fragmented transcripts.

## 4. Salmon SAF Quantification

- **Workflow**: Quantifies gene expression using Salmon's selective alignment (SAF), which is fast and accurate for transcript-level quantification.
- **Strengths**: Rapid quantification; effective for large datasets; less resource-intensive than alignment-based methods.
- **Limitations**: Relies on high-quality transcriptome annotation; may be less accurate for novel or poorly annotated genes.

## 5. Bowtie2 + RSEM Quantification

- **Workflow**: Uses Bowtie2 for alignment and RSEM for transcript quantification.
- **Strengths**: Accurate quantification; RSEM models read mapping uncertainty; suitable for both reference-guided and de novo approaches.
- **Limitations**: Slower than pseudo-alignment methods; requires careful parameter tuning for best results.

---

## Summary Table

| Method                    | Reference Needed    | Novel Transcript Discovery | Speed    | Accuracy | Typical Use Case           |
| ------------------------- | ------------------- | -------------------------- | -------- | -------- | -------------------------- |
| HISAT2 Ref Guided         | Yes                 | Yes (limited)              | Moderate | High     | Model organisms            |
| HISAT2 De Novo            | No                  | Yes                        | Moderate | Moderate | Non-model organisms        |
| Trinity De Novo           | No                  | Yes                        | Slow     | Moderate | Transcript discovery       |
| Salmon SAF Quantification | Yes (transcriptome) | No                         | Fast     | High     | Large-scale quantification |
| Bowtie2 + RSEM            | Yes/No              | Yes                        | Slow     | High     | Detailed quantification    |

---

## Recommendations

- Use **HISAT2 Reference Guided** for well-annotated genomes.
- Use **De Novo** methods (HISAT2 or Trinity) for non-model organisms or when reference is unavailable.
- Use **Salmon** for rapid quantification when transcriptome annotation is available.
- Use **Bowtie2 + RSEM** for high-accuracy quantification and uncertainty modeling.

---

For more details, see the respective method folders and scripts.
