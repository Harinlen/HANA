# Notes For Developer

## Current workflow

Now we expect to work in the following way:
- [x] `extract` to extract information from BAM or PAIRS files. (Input: `.fasta`, `.bam`/`.pairs`, Output: `.hmr_nodes`, `.hmr_reads`)
- [x] `draft` to generate graph information. (Input: `.hmr_nodes` and `.hmr_reads`, Output: nodes `.hmr_nodes` and edges `.hmr_edges`)
- [x] `partition` to generate group nodes. (Input: `.hmr_nodes` and `.hmr_edges`, Output: node group files `.hmr_group`)
- [x] `ordering` to sort the nodes. (Input: `.hmr_nodes`, `.hmr_edges` and `.hmr_group`, Output: node sequence files `.hmr_seq`)
- [x] `orientation` to decide the orientation of the node groups. (Input: `.hmr_nodes`, `.hmr_reads` and `.hmr_seq`, Output: chromosome sequence files `.hmr_chromo`)
- [ ] `build` to generate results. (Input: `.hmr_nodes`, `.fasta`, `.hmr_chromo`, Output: `.fasta`)