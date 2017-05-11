#!/bin/bash

set -eo pipefail

if [[ $TRAVIS_BRANCH != "master" && $TRAVIS_PULL_REQUEST == "false" && $TRAVIS_REPO_SLUG == "lcdb/lcdb-wf" ]]
then
    echo ""
    echo "Tests are skipped for non-master-branch pushes to the main repo."
    echo "If you have opened a pull request, please see the full tests for that PR."
    echo ""
    exit 0
fi

source activate lcdb-wf-test && snakemake -prs references.snakefile --configfile config/test_config.yaml --use-conda -j2
source activate lcdb-wf-test && snakemake -prs rnaseq.snakefile --configfile config/test_config.yaml --use-conda -j2
source activate lcdb-wf-test && py.test wrappers/test -n2 -v
