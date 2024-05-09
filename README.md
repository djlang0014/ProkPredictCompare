# APCCBSP24

Repo for Final Project.

NOTE! I forgot to use a virtual environment, so the requirements may install extraneous packages. My bad.

### About

This is a web tool that provides runs multiple prokaryotic gene
prediction tools, specifically Prodigal, MetaGeneAnnotator, FragGeneScan, and
Glimmer, then provides comparison metrics to evaluate their relative
performance/accuracy.

### Requirements

Requires valid installations of Prodigal, MetaGeneAnnotator,
FragGeneScan, and Glimmer. This can all be acquired from the tarball.

This is designed to run on Linux.

Python 3 is required. Python packages can be obtained by run "pip install -r requirements" after cloning this repo.

Final storage space is minimal, but there will be a
significant (~10mb) download during operation of new genome assemblies to
obtain GFF and FASTA files. These will be deleted after the results are
finalized.
