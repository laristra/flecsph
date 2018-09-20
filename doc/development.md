# General development guidelines

We discourage development on the `master` branch. 

Indeed, FleCSPH is an open-source project which is developed for a broad 
community of computer scientists, applied mathematicians and physicists
in various institutions across the globe.
Our ultimate goal is to provide a functional toolkit which is easy to 
download and deploy anytime, such that the code can be used and not broken 
due to various experimental features in the course of development.
That is why we use git version control and github.com platform as a tool 
for productive collaboration and as means to ensure code correctness and 
reliability.
Therefore, if you are directly working on the `master` branch, it can
interfere with the work of others.
Unless the changes are obvious typos or trivial quick fixes, all work and 
development should be done on separate personal branches.

We follow the standard git development workflow.
[Here](https://guides.github.com/introduction/flow/) is a good description 
of the git workflow that we adopt, specifically used in github service.
We further describe some specifics of the FleCSPH development workflow:

# Development Flow

## 1. Create your own branch, based on master branch
There are many ways to create a branch; here is what we recommend:
```{engine=sh}
   git fetch --all         # fetches all the updates from the origin without merging
   git checkout master     # switch to master branch
   git pull origin master  # update master if necessary
   git checkout -b <your_branch>  # branch off of master
```

### Branch naming conventions
We follow the FleCSI branch naming convention:

`<category>/<your username>/<description>`

The `<category>` tells what is the major objective of this branch. We are 
currently using `fix`, `feature`, `doc`, `stable` etc:
- _feature_: if you are adding a feature;
- _fix_: if you are working on a fix;
- _doc_: if you are editing documentation;
- _stable_: a stable, thoroughly tested version which can 
  be converted into a release.

With these naming conventions, it is easy to identify who is working on the
branch and what they are doing. When naming your new branch, take time to 
come up with a short and snappy description which effectively communicates 
the nature of your work. 

Also, please do not work on other people's branches unless you discussed it with
those people in advance. Modifying/deleting someone's branch without notification 
can be very dangerous in many ways.

## 2. Working and developing in your branch
Once you have your branch, you can develop it by editing the code and producing
one or more commits. While you are working, you can publish your new commits by 
pushing your branch upstream:
```{engine=sh}
   <editing>
   <commiting>
   ...
   git push origin <your_branch>
```

## 3. Pull request
As soon as you have finished developing your feature, fixing a problem, or otherwise
consider your work complete, you may want to request it to be added to the `master` 
branch. You can do it by creating a _pull request_ and assigning other developers to
examine your changes for approval.

Detailed introduction into the concept of pull requests can be found
[here](https://help.github.com/articles/creating-a-pull-request/). 
You can create pull requests both if you are an authorized LA Ristra developer,
and if you have created your own fork of the project.

When creating a pull request, you will need to assign reviewers and approvers. 
Please refrain from assigning yourself unless you have a split personality with 
one of your alter egos being a disengaged expert without conflict of interest.

## 4. Merge conflicts
Most of the time, your code development will generate merge conflicts. However,
if merge conflicts arise, it is better to contact the person who developed a 
conflicting path for correct conflict resolution.
If the conflicts are trivial, you can try to resolve them in your local git
repository and then push the merged version:
```{engine=sh}
   git checkout master
   git pull origin master # get recent updates
   git checkout <your_branch>
   git merge master # this pulls recent changes in master into your branch
   <... ! merge conflict ...>
   <resolve merge conflict by editing conflicted files>
   git add <resolved conflict files>
   git commit # complete merge after conflict has been resolved
   git push origin <your_branch>
```
[This link](https://help.github.com/articles/resolving-a-merge-conflict-using-the-command-line/)
gives further explanation about how to resolve merge conflicts.

Once you are done with your local branch and do not plan to use it anymore,
please delete it to avoid redundancies. You can simply do it via:
```{engine=sh}
   git branch -D <your_branch>
   git push --delete origin <your_branch> # or:
   git push origin :<your_branch>         # also deletes the branch on remote
```

# C++ development style guide

FleCSPH follows the FleCSI coding style, which in turn follows (in general) the Google coding conventions.
FleCSI coding style is documented here:
https://github.com/laristra/flecsi/blob/master/flecsi/style.md


# Contact

If you have any questions or concerns, please contact Oleg Korobkin (korobkin@lanl.gov), 
Julien Loiseau (julien.loiseau@univ-reims.fr), or Hyun Lim (hylim1988@gmail.com) 
