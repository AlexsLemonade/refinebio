# R dependencies
In this directory, we store all the scripts for installing our R dependencies.
These scripts are managed via `cranlock`, our system for locking versions of
transient dependencies. 

## Layout
In each directory, there are three files: `packages.txt`,
`versions.tsv`, and `dependencies.R`. `packages.txt` is the input file to
`cranlock`, and it holds the names of all of our R dependencies. `versions.tsv`
is part of the output of `cranlock`, and it holds all of R package versions for
the docker container that we are installing these packages in, and it is used
by `cranlock` to resolve the versions to install. `dependencies.R` is the script
that installs the dependencies in order as part of our docker build process.

## Usage
The first step to installing or adding new R dependencies is to install them
manually in a docker container to verify that all the versions are compatable.
This installation needs to take place inside the `Dockerfile` so that the
package versions can be queried by `cranlock`. Then, you can run `cranlock`,
passing in the path to `packages.txt` and the name of the docker container, and
it will output a `versions.tsv` file and a `dependencies.R` file in the same
directory as `packages.txt`.
