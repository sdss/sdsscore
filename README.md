# sdsscore_test

This repository contains files requires for SDSS-V operations at APO and LCO.

## Summary files

The directories `apo/summary_files` and `lco/summary_files` contain `confSummray` and `confSummaryF` files generates at the observatories for each configuration and FVC loop respectively. The structure is `<observatory>/summary_files/NNNXXX/NNNMXX/` where summary files are grouped in directories of one hundred and those grouped in directories that contain a thousand summaries. Each thousand grouping directory (`NNNXXX`) is a submodule pointing to a different repository in the form `sdss/sdsscore_apo_summary_files_NNNXXX` (e.g., [sdsscore_apo_summary_files_015XXX](http://github.com/sdss/sdsscore_apo_summary_files_015XXX)). This structure prevents the main sdsscore repository to grow larger than the limits allowed by GitHub.

### Cloning the `sdsscore_test` product

`sdsscore_test` can be cloned as any other git repository git

```console
git clone git@github.com/sdsscore_test
```

This will checkout the `main` branch and fetch any files that directly belong in the repository, but the summary files directories will appear empty. To initialise them do

```console
git submodule init
```

(this will take a long time as many GBs of data will need to be downloaded). Then you can update the submodules with

```console
git submodule update --recursive
```

The submodule directories will appear as detached repositories, pointing to the latest commits. If you want to commit files or make changes you will need to checkout the `main` branch in the submodule and make and commit changes there.

In general, refer to the [submodules](https://www.git-scm.com/book/en/v2/Git-Tools-Submodules) documentation.

### Adding a new submodule to `summary_files`

Every now and then the `confSummary` files being generated will catch up with the available `summary_files` submodules. If you need to add a new submodule, follow the following steps:

- Create a [new repository](https://github.com/new) in the SDSS organisation with name `sdsscore_OBS_summary_files_NNNXXX` (e.g., `sdsscore_apo_summary_files_015XXX`). Make sure you select the option to add a README file so that the resulting repository has a `main` branch. (We'll use `sdsscore_apo_summary_files_015XXX` as the example for the rest of these instructions).

- SSH to the observatory where you need to add the submodule. SSH as your own user and make sure you are forwarding an SSH key that also allows you to clone GitHub repositories.

- `cd /home/sdss5/software/sdsscore_test`. Ideally you'll be adding the submodule before the previous submodule is full (in that case ignore the rest of this point), but if that's not the case [jaeger](https://github.com/sdss/jaeger), the product that creates the `confSummary` files will have continued creating summary files and adding them to `sdsscore_test`. Some, but not all of the files may have been committed directly to `sdsscore_test` instead to a submodule. **Do not remove those files** as you may lose them and with them the associated observations. Instead move the directory you want to replace with a submodule to a location outside `sdsscore_test`: `mv apo/summary_files/015XXX ~/`. Stage and commit the changes, if any, to have a clean directory tree:

    ```console
    git add .
    git commit -m "Clearing apo/015XXX"
    git push
    ```

- Now you can add the new submodule with

    ```console
    git submodule add git@github.com:sdss/sdsscore_apo_summary_files_015XXX.git ./apo/summary_files/015XXX
    ```

- Run the following command to commit the changes. This is the same script that gets run by a cronjob hourly to commit new summary files

    ```console
    bash /home/sdss5/config/cronjobs/sdsscore.sh
    ```

    You will see a warning like this

    ```console
    warning: adding embedded git repository: apo/summary_files/015XXX
    hint: You've added another git repository inside your current repository.
    hint: Clones of the outer repository will not contain the contents of
    hint: the embedded repository and will not know how to obtain it.
    hint: If you meant to add a submodule, use:
    hint:
    hint:    git submodule add <url> apo/summary_files/015XXX
    hint:
    hint: If you added this path by mistake, you can remove it from the
    hint: index with:
    hint:
    hint:    git rm --cached apo/summary_files/015XXX
    hint:
    hint: See "git help submodule" for more information.
    [main be0a166d] Update submodules
    2 files changed, 4 insertions(+)
    create mode 160000 apo/summary_files/015XXX
    ```

    That is expected the first time. After the command finishes you should have a clean git session; running `git status` should return

    ```console
    On branch main
    Your branch is up to date with 'origin/main'.

    nothing to commit, working tree clean
    ```

    You can return the same `bash /home/sdss5/config/cronjobs/sdsscore.sh` and now it should complete without warnings.

- Now SSH to Utah and become the `sdssunit` user with `sudo su - sdssunit` (you need to belong to the group of users allowed to become `sdssunit`).

- Run `module load sdsscore/test` and `cd $SDSSCORE_DIR`.

- Run `git pull`, which should update `.gitmodules` and add the new submodule. Initialise the submodule with `git submodule init`. Now `cd` to the new submodule directory (`cd apo/summary_files/015XXX`) and run `git status`. If the submodule is not tracking the `main` branch (e.g. it's detached) do `git checkout main`.

- Repeat these steps to add another submodule.
