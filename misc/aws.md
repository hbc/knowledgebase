- no dirs in AWS, those strings are just prefixes
- use `--dryrun` to test a command

|Problem        | Unix rationale         | AWS spell                                                                     |
|---------------|------------------------|-------------------------------------------------------------------------------|
|copy files     |`cp * destination_dir`  |`aws s3 sync . s3://bucket/dir/`                                               |
|get file sizes |`ls -lh`                |`aws s3 ls --human-readable s3://bucket/dir/`                                  |
|copy bam files |`cp */*.bam /target_dir`|options order matters! <br/> `aws s3 sync s3://bucket/dir/ . --exclude "*" --include "*.bam"`|
