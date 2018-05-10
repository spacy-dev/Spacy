#!/bin/bash
set -e # Exit with nonzero exit code if anything fails

SOURCE_BRANCH="master"
TARGET_BRANCH="gh-pages"

# Pull requests and commits to other branches shouldn't try to deploy, just build to verify
if [ "$TRAVIS_PULL_REQUEST" != "false" -o "$TRAVIS_BRANCH" != "$SOURCE_BRANCH" ]; then
    echo "Skipping deployment of documentation."
    exit 0
fi

# Save some useful information
REPO=`git config remote.origin.url`
SSH_REPO=${REPO/https:\/\/github.com\//git@github.com:}
SHA=`git rev-parse --verify HEAD`

# Clone repository and checkout 'gh-pages'
git clone $REPO SpacyDoc
cd SpacyDoc
git checkout $TARGET_BRANCH
# Remove old documentation
rm -rf doc/*

# Add new documentation
mv  $TRAVIS_BUILD_DIR/build/doc/* doc/

# Set git user
git config user.name "Travis CI"
git config user.email "$COMMIT_AUTHOR_EMAIL"

# If there are no changes to the compiled out (e.g. this is a README update) then just bail.
if [ -z `git diff --exit-code` ]; then
    echo "No changes to the output on this push; exiting."
    exit 0
fi

# Add new documentation to git
git add --all doc/*
git commit -m"Update documentation for commit: ${SHA}"

chmod 600 travisci_rsa
eval `ssh-agent -s`
ssh-add travisci_rsa

# Push new documentation
git push ${REPO/https:\/\/github.com\//git@github.com:} $TARGET_BRANCH

