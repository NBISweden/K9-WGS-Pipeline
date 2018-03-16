#!/bin/bash

# Only auto merge develop branch for now
if [ "$TRAVIS_BRANCH" != "develop" ]; then
    exit 0;
fi

export GIT_COMMITER_EMAIL='johan.viklund@nbis.se'
export GIT_COMMITER_NAME='Johan Viklund'

git checkout master || exit
git merge "$TRAVIS_COMMIT" || exit
git push -q https://$GITHUB_TOKEN:x-oauth-basic@github.com/viklund/TestingAPIPushKey.git HEAD:master
