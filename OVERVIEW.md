# ImmunoScore: A Plain-Language Overview

*An in-depth explanation of what this project is, what it does, and why it matters, written so that anyone with a basic science background can follow it. No prior knowledge of cancer biology or programming is assumed.*

Kalixte Petrof, 2026.

---

## Table of contents

1. The thirty-second version
2. The problem: cancer and the immune system
3. Immunotherapy: releasing the brakes
4. The unsolved problem this project targets
5. The question I asked
6. A crash course in the biology you need
7. The data: single-cell RNA sequencing
8. What I actually did, step by step
9. The concepts behind the numbers
10. The results, told honestly
11. The biology of the finding
12. What could be wrong
13. Why it matters anyway
14. What I built
15. What comes next
16. Glossary
17. How to reproduce it yourself

---

## 1. The thirty-second version

Some cancer drugs work by waking up the immune system so it attacks the tumor. They are remarkable for a minority of patients and useless for the rest, and there is no reliable way to tell beforehand who is who. I used public data to build a "signature," a short list of genes whose activity in a certain immune cell tracks with whether a patient responds. I built it from melanoma patients, then tested it on a completely different cancer it had never seen. It still told responders apart from non-responders, though the sample was small enough that I treat the result as a promising lead rather than a proven fact. Everything is reproducible, and I measured exactly how uncertain the result is instead of hiding it.

---

## 2. The problem: cancer and the immune system

Your body is made of cells. Cells divide to grow and to repair, and normally they stop dividing when they should. Cancer begins when a cell breaks those rules. It keeps dividing, ignores the signals that tell it to stop, and eventually forms a tumor that can spread.

You already have a defense system against this: the immune system. Among its many cell types are T cells, which act like patrolling soldiers. A particular kind, the CD8 T cell, is a trained killer. When a CD8 T cell recognizes a cell as dangerous, infected by a virus or turned cancerous, it can destroy it directly. In principle, your immune system should be able to find and kill cancer cells.

Tumors survive anyway, because they cheat. One of their tricks is to press a molecular "off switch" on the T cells that come to attack them. The T cell arrives ready to kill, the tumor flips the switch, and the T cell stands down. The tumor grows in plain sight of an immune system that has been talked into ignoring it.

---

## 3. Immunotherapy: releasing the brakes

For most of medical history, cancer treatment meant attacking the tumor directly: cut it out, burn it with radiation, poison it with chemotherapy. Immunotherapy is different. Instead of attacking the tumor, it re-arms your own immune system.

The specific drugs this project is about are called checkpoint inhibitors, and the most common target is a switch named PD-1. PD-1 sits on the surface of T cells. It is one of the off switches tumors abuse. A drug called an anti-PD-1 antibody physically blocks that switch, so the tumor can no longer use it. The T cell stays active. It does the job it came to do.

When this works, it can work spectacularly. Some patients with advanced, previously untreatable cancers have lasting responses, occasionally years long. This is one of the most important advances in cancer medicine in decades, and it won a Nobel Prize in 2018.

Here is the catch, and it is the reason this project exists. Anti-PD-1 helps only a minority of patients, often somewhere between one in five and two in five depending on the cancer. The rest get the drug, get the side effects, get the bill, and lose months they did not have, with no benefit. And doctors usually cannot tell in advance which group a patient is in.

---

## 4. The unsolved problem this project targets

The question is simple to state and hard to answer: **before treatment starts, can you predict whether a given patient will respond to anti-PD-1?**

If you could, you would give the drug to the people it will help, and steer everyone else toward something with a better chance. You would save time, money, and suffering. Researchers around the world are working on this. It is an open problem, not a solved one.

I am a high-school student. I am not going to solve it. What I wanted to find out was how far someone with public data, free software, and persistence could get on a real version of this question, and whether I could do it honestly.

---

## 5. The question I asked

I narrowed the big problem to one specific, testable question:

> Can a pattern of gene activity in CD8 T cells, learned from patients with one cancer, predict treatment response in patients with a *different* cancer?

The reason this version is interesting: if a pattern only works in the exact cancer it was built on, it might just be memorizing that dataset. If the same pattern carries over to a different cancer, that hints at a shared biology of response, something more fundamental than one tumor type. That "carrying over" is called generalization, and it is the heart of the project.

---

## 6. A crash course in the biology you need

To follow what I did, you need four ideas. They build on each other.

**Cells and genes.** Almost every cell in your body carries the same complete set of instructions, your DNA, divided into about 20,000 genes. A gene is a recipe, usually for a protein, and proteins do most of the actual work in a cell.

**Gene expression.** A cell does not use all 20,000 genes at once. At any moment it switches some on and most off. A gene that is "on" gets copied into a short working message called mRNA, which the cell uses to build the corresponding protein. The set of genes a cell has switched on, and how strongly, is called its gene expression, and it defines what the cell is doing right now.

**Cell identity is expression.** A skin cell and a T cell carry the same DNA. They are different because they express different genes. So if you can read which genes a cell has switched on, you can tell what kind of cell it is and what state it is in.

**CD8 T cells.** These are the killer T cells, the ones anti-PD-1 is trying to unleash. They are the cells I focus on, because if any cells carry a signal about immunotherapy response, it is most likely these.

---

## 7. The data: single-cell RNA sequencing

To study gene expression you have to measure it, and the measurement technology matters a lot here.

**The old way: bulk sequencing.** You take a chunk of tumor, grind it up, and measure the average gene expression across all the cells mixed together. The problem is obvious once you see it: a tumor is a crowd of many cell types, and an average hides the individuals. It is like grading a class by the average of everyone's exam. You learn nothing about any single student.

**The new way: single-cell RNA sequencing (scRNA-seq).** This technology separates the tissue into individual cells and reads the gene expression of each one separately. Now instead of the class average you have every student's exam. You can see the killer T cells on their own, apart from the tumor cells, the support cells, and everything else.

**What the data physically is.** The result is an enormous table. Each row is one cell, each column is one gene, and each number is how active that gene was in that cell. A single dataset can have tens of thousands of cells and tens of thousands of genes. Most of the table is zeros, because any given cell only expresses a fraction of all genes at a detectable level. Computer scientists call this a sparse matrix, and there are efficient ways to store and work with it, which the code uses.

**The two datasets I used.** Both are public, published by professional labs, and both track real patients through anti-PD-1 treatment with a recorded outcome (responder or non-responder).

- *Melanoma* (skin cancer), from Sade-Feldman and colleagues, 2018. About 16,000 immune cells across roughly 48 samples. I used this to build the signature.
- *Basal cell carcinoma* (a different skin cancer), from Yost and colleagues, 2019. About 53,000 cells across 11 patients with response labels. I used this only to test the signature, never to build it.

Keeping those two roles strictly separate, build on one, test on the other, is what makes the test honest.

---

## 8. What I actually did, step by step

The whole project is a pipeline: data goes in one end, a tested result comes out the other. Here is every stage in plain language, and why each one is there.

**Step 1: Load and quality-control the cells.** Raw single-cell data is full of junk: droplets that caught no real cell, and cells that were dying when measured. Dying cells leak their contents and show a tell-tale signature of high mitochondrial gene activity. I filter those out. If you skip this, noise pollutes everything downstream.

**Step 2: Normalize.** Some cells are measured more deeply than others, just by chance of the technology. Raw counts are therefore not comparable cell to cell. Normalizing rescales every cell to a common total, then applies a logarithm to keep a handful of very active genes from drowning out the rest. After this, a number means the same thing in every cell.

**Step 3 (a practice run): cluster and identify cell types.** Before touching the real question, I ran the standard workflow on a tiny built-in dataset and let the computer group cells by similarity. The groups it found matched known immune cell types exactly, marked by the genes immunologists already use (for example CD3D for T cells, CD79A for B cells). That match was the proof my pipeline was working correctly, before I trusted it on anything important.

**Step 4: Focus on baseline CD8 T cells.** From the melanoma data I kept only the CD8 T cells, and only from samples taken *before* treatment started. "Before treatment" matters: a predictor you can only compute after the drug is useless. I want a signal present at baseline, when the decision is actually made.

**Step 5: Find the difference.** Among those baseline CD8 T cells, I compared the responders to the non-responders and asked which genes were more active in the people who went on to respond. This comparison is called differential expression, and it produces a ranked list of candidate genes.

**Step 6: Clean the list.** The top of that raw list was polluted with genes that are famous technical artifacts in this kind of data: ribosomal genes, mitochondrial genes, and stress-response genes that switch on simply because the cells were handled. They carry no real biological signal about response. I wrote a filter to drop them, then took the top genes that remained.

**Step 7: The signature.** What survives is the signature: fifty genes whose combined activity is associated with response. Think of it as a fingerprint of a "responsive" immune state.

**Step 8: The real test, on a different cancer.** I took that melanoma signature and, without changing it at all, used it to score every CD8 T cell in the basal cell carcinoma data. Scoring just asks, for each cell, how strongly it expresses the signature genes as a group. I then averaged each patient's cells into one number per patient, and asked whether responders scored higher than non-responders. This is where the AUC of 0.77 comes from, explained in the next section.

**Step 9: A more rigorous version (pseudobulk).** There is a statistical trap here. Cells from the same patient are not independent measurements; they are the same person sampled many times. Treating 2,000 cells from 19 people as 2,000 independent data points overstates your evidence. The fix is to collapse each patient's cells into a single averaged profile first, so the unit of analysis is the patient, not the cell. I rebuilt the signature this way. Interestingly, my first attempt at it scored *worse* than random, because at the patient level a few barely-expressed genes produced flukishly strong statistics. I diagnosed that, added a filter requiring genes to be genuinely expressed, and the rigorous version then scored slightly better than the original. That debugging story, getting it wrong, understanding why, and fixing it, is part of the work.

**Step 10: Measure the uncertainty.** A single number like 0.77 means little without knowing how shaky it is. I ran two standard tests. A bootstrap resamples the patients thousands of times to produce a confidence interval, the range the true value probably lives in. A permutation test shuffles the response labels thousands of times to see how often pure chance produces a score this good. Both are explained below. Both told the same story: real signal, but uncertain, because eleven patients is not many.

---

## 9. The concepts behind the numbers

**A gene signature** is a list of genes treated as a group. Instead of asking about one gene, you ask whether a cell or a patient has the whole set turned up together.

**Scoring** turns that list into a single number per cell: roughly, how far above the background this cell's signature genes are, taken together. Higher score means the cell looks more like the "responsive" state.

**AUC (area under the ROC curve)** measures how well the score separates two groups, on a scale from 0.5 to 1.0. The cleanest way to read it: pick one responder and one non-responder at random; the AUC is the probability that the responder has the higher score. 0.5 means you are guessing (a coin flip). 1.0 means perfect separation every time. My result, 0.77, means about a 77 percent chance the responder ranks higher. Good, clearly better than chance, not perfect.

**A p-value** answers a specific question: if the signature actually meant nothing, how often would pure luck produce a result at least this good? Mine is 0.08, which is 8 percent. By long-standing convention, scientists call a result "statistically significant" only when that probability drops below 5 percent. Mine is above that line, so the honest word is "suggestive," not "significant." It could still be a fluke, and I say so.

**A confidence interval** is the range the true value most likely falls in. Mine runs from 0.40 to 1.00. That is a wide range, and it dips below 0.5, which is another way of saying that with only eleven patients I cannot rule out that the real performance is no better than chance. The width is the honesty.

**A bootstrap** is how you get that interval when you do not have a formula handy. You pretend your eleven patients are the whole world, draw eleven of them at random with replacement, recompute the score, and repeat a few thousand times. The spread of those repeats estimates how much your answer would wobble if you had sampled different patients.

**A permutation test** is how you get the p-value. You scramble the response labels so they are random, recompute the score, and repeat thousands of times to build a picture of what pure chance looks like. Then you see where your real score falls in that picture. Mine sits at the right edge of the chance distribution, close to it, not far beyond it.

**Pseudoreplication and pseudobulk.** Pseudoreplication is the trap from Step 9: counting non-independent measurements as if they were independent. Pseudobulk is the fix: average each patient down to one profile so each patient counts once.

---

## 10. The results, told honestly

Here is the whole result in one place.

- The melanoma-built signature, scored on basal cell carcinoma, separated responders from non-responders at a patient-level AUC of about 0.77.
- That estimate carries a wide 95 percent confidence interval, 0.40 to 1.00, and a permutation p-value of 0.08.
- So: a real, encouraging signal, but not statistically significant at the usual threshold, because the validation set had only eleven patients.

I want to be precise about what is solid and what is not. The *prediction* is suggestive and needs more data. The *biology* underneath it is solid, because the signature is built from genes that immunologists already associate with response, and I recovered them from raw data without being told to. Recovering a known result from scratch is strong evidence the method works, even when the headline number is uncertain.

A weaker project would report 0.77 and stop. Reporting the confidence interval and the p-value next to it, and calling the result a lead rather than a discovery, is the part I am proudest of.

---

## 11. The biology of the finding

The signature is led by a recognizable group of genes: TCF7, CCR7, IL7R, CD28, SELL. These mark a particular state of CD8 T cell, often called stem-like or memory-like.

The intuition: not all killer T cells are the same. Some are exhausted, worn out from constant fighting, near the end of their useful life. Others are stem-like, fresh, long-lived, able to renew themselves and produce new fighters. A tumor where the immune response still has these fresh, self-renewing cells at baseline has reinforcements available when the drug releases the brakes. A tumor whose T cells are already exhausted has no reserves to call on. Patients with more of the stem-like cells tend to respond, and this matches published findings. The drug works best where there is still an army left to un-pause.

The more rigorous pseudobulk version of the signature highlighted a slightly different but related group, genes involved in T-cell activation and signaling (LAT, DGKA, PTK2B, CD28, CD69). Both stories point the same way: a functional, ready immune state at baseline predicts response. That two different analyses landed on the same theme is reassuring.

---

## 12. What could be wrong

Good science states its own weaknesses first. Here are mine, in order of seriousness.

**Small sample.** Eleven validation patients is the headline limitation. It is why the confidence interval is wide and the p-value is above the significance line. Everything else is secondary to this.

**Batch effects.** The two datasets were produced by different labs, at different times, with different equipment. Some of the difference between them is technical, not biological. It is possible that part of what looks like a transferable signal is actually a technical fingerprint of one dataset. I have not fully ruled this out, and saying so is part of being honest.

**The signature may be unspecific.** It is possible the genes track general T-cell health rather than something specific to anti-PD-1 response. A healthier immune system might look better on many measures. Distinguishing "specifically predicts this drug" from "generally looks healthy" would need more careful controls.

**It is a reanalysis.** I did not generate new patient data. I reanalyzed datasets other scientists collected and published. That is legitimate and common, but it means the project's ceiling is set by what those datasets contain.

**Overfitting risk.** Any time you search thousands of genes for the ones that separate two groups, you risk finding patterns that are specific to your particular patients rather than general truths. The cross-cancer test partly guards against this, which is exactly why I built the project around it, but small samples make the risk real.

**It is not medical advice.** Nothing here is validated for clinical use. It does not guide treatment for any real patient. It is a research and learning project.

---

## 13. Why it matters anyway

Given all those caveats, why is this worth anything?

Because the value is not a clinical claim. The value is threefold. First, it shows that a serious, modern question in cancer immunology can be approached with public data and open tools by someone willing to learn, and that the method genuinely recovers real biology. Second, the specific idea, that a response signature might transfer across cancer types, is a real and testable hypothesis, and this project is a small piece of evidence for it. Third, the whole thing is reproducible and honestly reported, which is the part of science that actually matters and the part most often skipped.

In other words: it is a real, if small, contribution, executed and reported the way real science is supposed to be.

---

## 14. What I built

The project is not just an analysis. It is a set of working pieces:

- A reproducible **pipeline** of short scripts that take raw single-cell data and produce the signature, the validation, the figures, and the statistics.
- The **signature files** themselves, the actual gene lists, saved so anyone can inspect or reuse them.
- **Unit tests** that check the code's core logic behaves correctly, so the work can be trusted and rerun.
- A small **interactive tool** that scores single-cell data against the signature.
- A **website** that explains the finding for a general reader.
- A public, **open-source repository** with everything needed to reproduce every number.

---

## 15. What comes next

The clearest next step follows directly from the main limitation. The result is uncertain because eleven patients is too few, so the fix is more patients. Several more public anti-PD-1 single-cell datasets exist for other cancers. Pooling them would increase the sample size, tighten the confidence interval, and show whether the signal holds up or fades. After that: a written report submitted to a peer-reviewed venue, feedback from a researcher in the field, and a science competition where experts evaluate the work.

---

## 16. Glossary

- **Anti-PD-1 / checkpoint inhibitor:** a drug that blocks an immune "off switch" so T cells stay active against the tumor.
- **AUC:** a 0.5-to-1.0 score for how well a signature separates two groups; 0.5 is chance, 1.0 is perfect.
- **Batch effect:** a technical difference between datasets that can be mistaken for a real biological difference.
- **Bootstrap:** repeatedly resampling your data to estimate how uncertain a result is.
- **CD8 T cell:** the killer immune cell that anti-PD-1 aims to unleash.
- **Confidence interval:** the range a true value most likely falls in.
- **Differential expression:** finding genes that are more active in one group than another.
- **Gene expression:** which genes a cell has switched on, and how strongly.
- **Pseudobulk:** averaging each patient's cells into one profile so each patient counts once.
- **Permutation test:** shuffling labels many times to see what pure chance looks like, giving a p-value.
- **p-value:** the probability of seeing a result this good by luck if the effect were not real.
- **Signature:** a list of genes treated as a group whose combined activity tracks with something.
- **Single-cell RNA sequencing:** technology that reads gene expression in each cell separately.
- **Stem-like / memory T cell:** a fresh, long-lived, self-renewing T cell state linked to immunotherapy response.

---

## 17. How to reproduce it yourself

Everything is in the repository. With Python installed, you create an environment, install the listed packages, download the two public datasets from TISCH2 (the README gives the exact names and where to put them), and run the scripts in order. Each script prints what it is doing. The same numbers in this document come out the other end. That is the point of building it this way: you do not have to take my word for any of it.
