# cWB pipeline public libraries

This repository contains all source files to
to compile cWB pipeline public libraries 

## Cloning the cWB library public repository

    git clone git@gitlab.com:gwburst/public/library.git

## Adding a public branch (called "public") to an existing repository

### Adding a new remote (called "public") and fetching all branches
    git remote add public git@gitlab.com:gwburst/public/library.git
    git fetch --all
### Creating locally  a new branch called "public" taken from the "public" remote     
    git checkout -b public public/public
