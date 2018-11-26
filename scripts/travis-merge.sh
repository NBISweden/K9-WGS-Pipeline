#!/bin/bash
set -e

# Only auto merge develop branch
if [ "$TRAVIS_BRANCH" != "develop" -o "$TRAVIS_PULL_REQUEST" != "false" ]; then
    exit 0;
fi
exit 0

git config user.email "johan.viklund@nbis.se"
git config user.name "Travis CI"

git fetch --depth=50 origin refs/heads/master:master
git checkout master
git merge --ff-only "$TRAVIS_COMMIT"
git push -q https://$GITHUB_TOKEN:x-oauth-basic@github.com/NBISweden/K9-WGS-Pipeline.git master
