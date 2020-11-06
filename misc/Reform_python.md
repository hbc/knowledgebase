# Using Reform to create custom genome
> https://gencore.bio.nyu.edu/reform/

- You can find good example with cellRanger in this link.


## Usage guide
### Conda env. and install
```
module load gcc/6.2.0 conda2/4.2.13
conda activate python_3.6.5 # you need an environment with python3 activated. (this one is custom for me)
pip3 install biopython
git clone https://github.com/gencorefacility/reform.git 
cd reform/
```

### The command code and the files we need.
```
  --chrom=<chrom> \
  --position=<pos> \ 
  --in_fasta=<in_fasta> \
  --in_gff=<in_gff> \
  --ref_fasta=<ref_fasta> \
  --ref_gff=<ref_gff>
```
- **chrom** ID of the chromsome to modify

- **position** Position in chromosome at which to insert <in_fasta>. Can use -1 to add to end of chromosome. Note: Either position, or upstream AND downstream sequence must be provided.

- **upstream_fasta** Path to Fasta file with upstream sequence. Note: Either *position*, or *upstream AND downstream* sequence must be provided.

- **downstream_fasta** Path to Fasta file with downstream sequence. Note: Either *position*, or *upstream AND downstream* sequence must be provided.

- **in_fasta** Path to new sequence to be inserted into reference genome in fasta format.

- **in_gff** Path to GFF file describing new fasta sequence to be inserted.

- **ref_fasta** Path to reference fasta file.

- **ref_gff** Path to reference gff file.

Example:
```ruby
python3 reform.py 
 --chrom=X \
 --position=3 \
 --in_fasta=in.fa \
 --in_gff=in.gff3 \
 --ref_fasta=ref.fa \
 --ref_gff=ref.gff3
 ```
 
> in.fa
 ```
 >input_sequence
TGGAGGATCG
```

> ref.gff
<img src="https://github.com/yoonsquared/knowledgebase_image/blob/master/gff3.png">

> in.gff
<img src="https://github.com/yoonsquared/knowledgebase_image/blob/master/in_gff3.png">

> reformed.gff
<img src="https://github.com/yoonsquared/knowledgebase_image/blob/master/reformed_gff.png">
