# Ska Ssc Low Scicomm Aa05

## Description

This repository is for software projects within the SKA Observatory that is not targeting construction. It is aimed for those who are not directly software developers

We try to keep the required steps a new non-developer user needs to take to a minimum so that you can start contributing as quickly as possible, but we will also add further information for those who are interested in following our best practices.

## Table of Contents

* [Required steps](#required-steps)
   1. [Naming your projects](#naming-your-projects)
   2. [Project license](#project-license)
   3. [Git basics](#git-basics)
* [Optional steps](#optional-steps)
    1. [Git branch based development](#git-branch-based-development)
    2. [Tracking branches and commit with Jira](#tracking-branches-and-commit-with-jira)
    3. [Continuous Integration and Deployment](#continuous-integration-and-deployment)

## Required steps

### Naming your projects

Projects must follow the follow the SKAO naming convention defined on Software Architecture Decision [ADR-25](https://confluence.skatelescope.org/display/SWSI/ADR-25+General+software+naming+convention).

Repository names shall clearly map to a particular element of the SKA software architecture, as described in the SKA software design documentation. That is to say, someone familiar with the SKA software architecture should be able to identify the content of a repository just by its name.

Names shall be all lowercase starting with “ska-” whereas multiple words shall be separated by hyphens.

### Project license

SKA organization promotes a model of open and transparent collaboration. In this model collaboration is made possible using permissive licenses, and not by pursuing the ownership of code by the organization. Copyright will thus belong to the institutions contributing to source code development for the lifetime of the project, and software developed for the SKA telescope will be available to the wider community as source code. Redistribution of the SKA software will always maintain the original copyright information, acknowledging the original software contributor.

The repository will be setup with the SKAO BSD 3-clause LICENSE. Any exception to this shall be justified and agreed with SKA office Software Quality Assurance engineer.

### Git basics

If you are new to Git and GitLab, we recommend you to read the [Learn Git page on GitLab](https://docs.gitlab.com/ee/gitlab-basics/start-using-git.html). There you will find a comprehensive guide to Git, including how to create a repository, how to create branches, how to make commits, and how to open pull requests.

Still if you are in a hurry, for a quick start, you can follow the steps below.

#### Tell your collaborators who you are
This will set your name and email address in your local Git configuration.

```shell
git config --global user.name "John Doe"
git config --global user.email johndoe@example.com
```

#### Clone the repository
This will create a copy of the repository on your local machine.

```shell
git clone PROJECT_GITLAB_URL
```

#### Make changes
Make changes to the files in your project as required.

#### Stage your changes
This will add the changes you made to the staging area.

```shell
git add FILE_NAME
```

#### Commit your changes
This will save the changes you made to the local repository.

```shell
git commit -m "A short message describing the changes you made"
```

#### Push your changes
This will upload the changes you made to the remote repository.

```shell
git push
```

#### Update your local working copy
If the main branch has been updated, you will need to update your local working copy to include the latest changes done by your colleagues.

```shell
git pull
```

## Optional steps

This section presents optional functionality that can be included in your project.

### Git branch based development

Branch based development is a common practice in software development. It allows you to work on new features or bug fixes without affecting the main branch. This way you can experiment with new ideas and make changes without worrying about breaking the main branch.

Besides what was already presented in the previous section, you can also follow the steps below to create a new branch, publish your branch and open a merge request.

#### Create a new branch
This will create a new branch based on the current branch you are on.

```shell
git checkout -b BRANCH_NAME
```

#### Publish your branch
This will create the branch on the remote repository so that others can see it.

```shell
git push origin BRANCH_NAME --set-upstream
```

#### Push your changes to the branch
This will upload the changes you made to the remote repository.

```shell
git push origin BRANCH_NAME
```

#### Open a merge request
This will open a merge request so that your changes can be reviewed and merged into the main branch.

1. Go to the project's GitLab page.
2. Click on the "Merge Requests" tab.
3. Click on the "New merge request" button.
4. Select the source branch (the branch you created) and the target branch (the main branch).
5. Click on the "Compare branches and continue" button.
6. Add a title and a description for your merge request.
7. Click on the "Submit merge request" button.
8. Wait for your colleagues to review your changes.
9. Make any necessary changes based on the feedback you receive.
10. Once your changes have been approved, click on the "Merge" button.
12. Your changes have been merged into the main branch.
13. You can now delete the branch you created.

#### Update your branch
If the either your branch or the main branch has been updated, you will need to update your local working copy to include the latest changes done by your colleagues.

```shell
git checkout main
git pull
git checkout BRANCH_NAME
git merge main
```

### Tracking branches and commits with Jira

If you are using Jira to track your project, you can include the Jira issue number in your commit messages. This will allow you to link your commits to the Jira issue and track the progress of your project.

In order to use this integration you need the following:

1. Create a Jira ticket: e.g. `MT-123`
2. Create a branch with the reference from the jira ticket as the branch name in lowercase followed by a description of your work using dashes as spaces: e.g. `mt-123-implement-x-feature`
3. Make your changes and commit them with the Jira ticket reference in the commit message: e.g. `git commit -m "MT-123: Add missing file"`
4. Open a merge request and include the Jira ticket reference in the title: e.g. `MT-123: Implement X feature`

Now your commits and merge requests will be linked to the Jira ticket and you will be able to track the progress of your project.

### Continuous Integration and Deployment

If you want to include continuous integration in your project, you can use GitLab CI/CD. This will allow you to automatically build and test your project every time you make a change.

Please consult the [SKAO Developer Portal on CI/CD](https://developer.skao.int/en/latest/tools/ci-cd.html) for more information on how to set it up for your project.