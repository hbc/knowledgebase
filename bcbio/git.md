# Git tips

- [Pro git book](https://git-scm.com/book/en/v2)
- https://github.com/Kunena/Kunena-Forum/wiki/Create-a-new-branch-with-git-and-manage-branches
- https://nvie.com/posts/a-successful-git-branching-model/
- http://sandofsky.com/blog/git-workflow.html
- https://blog.izs.me/2012/12/git-rebase

# Sync with upstream/master, delete all commits in origin/master branch
```
git checkout master
git reset --hard upstream/master
git push --force
```

# Sync with upstream/master
```
git fetch upstream
git checkout master
git merge upstream/master
```

# Feature workflow
```
git checkout -b feature_branch
# 1 .. N
git add -A .
git commit -m "sync"
git push?

git checkout master
git merge --squash private_feature_branch
git commit -v
git push
# pull request to upstream
# code review
# request merged
git branch -d feature_branch
git push origin :feature_branch
```
