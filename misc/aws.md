Can be run from ` /n/app/bcbio/dev/anaconda/bin/aws` or from  `/usr/bin/aws` on the O2 transfer nodes.

## Setup S3 bucket
Setup your S3 bucket with: 
`foo/aws configure`
You will need your 
- AWS Access Key ID
- AWS Secret Access Key
- Default region name
- Default output format

    
## Interact with AWS bucket      

- no dirs in AWS, those strings are just prefixes
- use `--dryrun` to test a command

|Problem        | Unix rationale         | AWS spell                                                                     |
|---------------|------------------------|-------------------------------------------------------------------------------|
|copy files     |`cp * destination_dir`  |`aws s3 sync . s3://bucket/dir/`                                               |
|get file sizes |`ls -lh`                |`aws s3 ls --human-readable s3://bucket/dir/`                                  |
|copy bam files |`cp */*.bam /target_dir`|options order matters! <br/> `aws s3 sync s3://bucket/dir/ . --exclude "*" --include "*.bam"`|
