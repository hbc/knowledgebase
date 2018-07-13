# How to remove all files except the ones you want: 

First, expand rm capabilities by:
`shopt -s extglob`

Then use find and remove:
`find . ! -name 'file.txt' -type f -exec rm -f {} +`