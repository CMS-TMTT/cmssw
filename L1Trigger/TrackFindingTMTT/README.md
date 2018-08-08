# Setup instructions to just run

To checkout and run, or if your modifications won't need to be put into the central repository, do:

```
cmsrel CMSSW_9_3_8
cd CMSSW_9_3_8/src
cmsenv

git cms-init
git remote add -t portTMTT TMTT git@github.com:EmyrClement/cmssw.git
git fetch TMTT portTMTT
git cms-checkout-topic EmyrClement:portTMTT

scramv1 b -j 8

cd TMTrackTrigger/TMTrackFinder/test/
cmsRun tmtt_tf_analysis_cfg.py inputMC=../../MCsamples/937/RelVal/TTbar/PU200.txt Events=10
```

Note that if you want to checkout other packages e.g. the tracklet software, which contains the ntuple makes, you should replace ```git cms-checkout-topic``` with ```git cms-merge-topic```

# Setup instructions for making modifications

Follow the above instructions.  At this point, you should be on a local branch called portTMTT.

Below is a simple example of making modifications, pushing them to your remote repository, and making a pull request back to the repository you want your changes to end up in.  In this example, we will use the following as the "central" repository (something equivalent of trunkSimpleCode9 in svn) : https://github.com/EmyrClement/cmssw  Follow the link, and fork the repository to your own account.  You then need to add your newly forked repository as a remote repository in your local working area, which we will call ```origin```:
```
git remote add origin <url>
```
You can get the url by clicking on "Clone or download" on the webpage for YOUR repository, and will be something like ```git@github.com:<GitHubUsername>/cmssw.git```

Lets change branch to one called "myChanges":
```
git checkout -b myChanges
```
Modify some files:
```
echo "#Hello World" >> TMTrackTrigger/TMTrackFinder/test/tmtt_tf_analysis_cfg.py
```
Check the status of your modifications:
```
git status
```
Which should show something like:
```
# On branch myChanges
# Changes not staged for commit:
#   (use "git add <file>..." to update what will be committed)
#   (use "git checkout -- <file>..." to discard changes in working directory)
#
#	modified:   TMTrackTrigger/TMTrackFinder/test/tmtt_tf_analysis_cfg.py
#
no changes added to commit (use "git add" and/or "git commit -a")
```
You can undo (revert) your changes (as explained in the message for ```git status```) with ```git checkout -- <file1> <file2> ...```
To see your modifications, you can do:
```
git diff TMTrackTrigger/TMTrackFinder/test/tmtt_tf_analysis_cfg.py
#or just
#git diff
```
Which should show something like:
```diff --git a/TMTrackTrigger/TMTrackFinder/test/tmtt_tf_analysis_cfg.py b/TMTrackTrigger/TMTrackFinder/test/tmtt_tf_analysis_cfg.py
index 7b693b9..fc59e07 100644
--- a/TMTrackTrigger/TMTrackFinder/test/tmtt_tf_analysis_cfg.py
+++ b/TMTrackTrigger/TMTrackFinder/test/tmtt_tf_analysis_cfg.py
@@ -189,3 +189,4 @@ if options.outputDataset == 1:
   process.writeDataset.outputCommands.append('keep  *_TTAssociator*_TML1Tracks*_*')
 
   process.schedule = cms.Schedule(process.p, process.pa, process.pd)
+#Hello World
```
Add the files, and commit:
```
git add TMTrackTrigger/TMTrackFinder/test/tmtt_tf_analysis_cfg.py
git commit -m "Modification to tmtt_tf_analysis_cfg comfig file"
```
Doing ```git status``` will now show:
```
# On branch myChanges
nothing to commit, working directory clean
```
Now push these changes to your remote repository.
```
git push origin myChanges 
```
You are now ready to make a pull request back to the central repository.  Go the webpage for your remote repository, and you should see a box stating you have just pushed some changes to the myChanges branch, and gives the option to "Compare & pull request".  Make the pull request to merge your changes in the myChanges branch to EmyrClement/cmssw:portTMTT.  You shoulw end up with a pull request that looks like this : https://github.com/EmyrClement/cmssw/pull/1

# Pull changes from the central TMTT repository
If changes (new commits) have been made to the central repository since you first checked it out, you can rebase your branch.  If you have any local changes, you need to ```git add``` and ```git commit``` them first.  Then do:
```
git pull --rebase TMTT portTMTT
```
You may have to resolve conflicts, in the usual git way.

Alternatively, if you don't want to commit your local changes, you can do ```git stash```, do the rebase, and then ```git stash pop```.

# To run in a newer CMSSW release
The above repository was setup in CMSSW_9_3_8.  If you want to run our TMTT software in a newer release, e.g. CMSSW_10_2_0, the command above will attempt things like merging all the differences in all packages between CMSSW_9_3_8 and CMSSW_10_2_0.  This can be avoided by manually performing a sparse checkout and tell git to only consider the TMTT directory when performing a checkout (which is similar in part is what the above commands do, as you don't need to checkout all of CMSSW).

```
cmsrel CMSSW_10_2_0
cd CMSSW_10_2_0/src
cmsenv

git cms-init
git remote add -t portTMTT TMTT git@github.com:EmyrClement/cmssw.git
git fetch TMTT portTMTT

echo "/TMTrackTrigger/TMTrackFinder" >> .git/info/sparse-checkout
echo "/TMTrackTrigger/MCsamples" >> .git/info/sparse-checkout

git checkout -b myChanges TMTT/portTMTT
```
You can then make changes as before, and create a pull request back to the same repository as before.

# To make a new branch of our software for a new CMSSW release
If you want to make a new branch of our software for a new CMSSW release e.g. when we move permanently to a new CMSSW release, you can do the following.  This example is moving from CMSSW_9_3_8 (the release the portTMTT branch was created in) to CMSSW_10_2_0
```
cmsrel CMSSW_10_2_0
cd CMSSW_10_2_0/src
cmsenv

git cms-init
git cms-rebase-topic -o CMSSW_9_3_8 EmyrClement:portTMTT
```
You will end up on a branch called portTMTT, which will have our TMTT software on top of CMSSW_10_2_0.  You can change to a new branch e.g. portTMTT_10_2_0, and push the new branch to the remote repository (if you have permissions).  If this new CMSSW release and branch should become the "master" branch everyone should work from, then you should also update the recipes/documentation.

