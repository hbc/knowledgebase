# AWS

Remember there is no dirs in AWS, those strings are just prefixes.

|What you wanted to achieve| Your Unix intuition | How it is in AWS|
|-------------------------|---------------------|-----------------|
|copy all files from the current dir to aws|`cp * destination_dir` |`aws s3 sync . s3://yourbucket/yourdir/`|
  
