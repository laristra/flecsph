# General Logistic

FleCSPH is open-source code which is developed among many different
people from domain scientists to computer scientists in many different
institutions. We would like to use version control to collaborate
productively and ensure correct code. We also aim to adapt our toolkit 
for the community. Thus, stable `master` branch is required to provide 
the best performance for users.

Furthermore, if you are directly working on the `master` branch, it can
interfere with other people's work and it makes hard to maintain quality
of the FleCSPH.
Unless the chages are obvious typos or simple quick fixed, all major
works/developments should be done at different branch.

[Here](https://www.atlassian.com/git/tutorials/comparing-workflows) is
good description on git workflows.

So, we generally ask people who develop/maintain the FleCSPH to follow
below flow.

# Development Flow

If not otherwise indicated, the FleCSPH development style follows below. 

## 1. Create your own branch for development
There are many ways to create the branch but here is simple way:
```{engine=sh}
   git fetch --all # This fetchs all the updates from the origin without merging
   git checkout master
   git pull origin master
   git checkout -b <your_branch>
```

## 2. Working and developing in your branch
Once you create your branch, you can edit/commit as usual git enviroment such that:
```{engine=sh}
   <editing>
   <commiting>
   ...
   git push origin <your_branch>
```
and so on

## 3. Pull request
Once you introduce major changes and/or restructures of code, you
may want to merge it to `master`branch. In git verseion control
system, you can simply use `git merge` commands but we ask to use
`pull request`.

A detail instruction for pull request can be found
[here](https://help.github.com/articles/creating-a-pull-request/). Note
that if you are added in the project (i.e. authorized developer for LARISTRA),
please do not sign yourself for pull request.

Most of time, your own development will not be problematic when you
try to merge. Sometimes you may encounter with merging problem. If
that case is happend, please ask below contact person before you
request merge.

Once you are done with your local branch and no plan to use anymore,
please delete to avoid any possible redundancies. You can simply do it via:
```{engine=sh}
   git branch -D <your_branch>
   git push --delete origin <your_branch>
```

# Contact

If you have any questions or concern related with this, please contact Oleg Korobkin (korobkin@lanl.gov), Julien Loiseau (julien.loiseau@univ-reims.fr), or Hyun Lim (hylim1988@gmail.com) 
