Shared Python code for the Atmospheric Chemistry Research Group (University of Bristol)

***First Installation***
Click on the link in the automated email sent by Matt and set up your bitbucket username and password.

Then check out the code using your bitbucket username. By doing the below:

From the directory where you want the "acrg" directory to live (e.g. for me on my laptop it's /Users/as13988/Documents/Work/Python/), type: 

git clone https://YOUR_USERNAME@bitbucket.org/mrghg/acrg.git

When you check out the code bitbucket automatically creates the "acrg" directory and fails if there's already an "acrg" directory. So if you want to keep your old google based "acrg" directory rename it.


***Adding a file to the repo***
From your directory, first add a file to the SVN library:

git add YOUR_FILE_NAME

Then "commit" all the changes you've made:

git commit -m "ADD_SOME_DESCRIPTIVE_COMMENTS_HERE_IN_QUOTES"

Note that it won't work if you don't add the -m with some comments.

Then "push" it to the central repository (the "origin master" isn't really necessary, but apparently it's a good idea for when we start using branches).

git push origin master

***Update regularly***
To make sure you're up-to-date with everyone else's changes, regularly type (or set a cron job), from inside the folder:

git pull

***Changing someone else's code***
First, make sure you have the most up-to-date version of the repo (see above). Then, edit the code on your machine. Then, when you're done, just commit the changes back to the repo:

git add FILENAME

git commit -m "ADD_SOME_DESCRIPTIVE_COMMENTS_HERE_IN_QUOTES"

git push origin master

*Adding an ssh key*

Git will ask you for a password every time you want to connect using https. You can get around this by adding your RSA key to your account. Instructions are here: https://confluence.atlassian.com/display/BITBUCKET/Set+up+SSH+for+Git. You can probably start at step 6.

***Conflicts***

Conflicts are inevitable. Git sorts this out for you. If you try to "push" changes from your local repository when you're not up-to-date, it'll complain. In that case, pull the latest changes from the server. There are now a few things that can happen:

1. The "pulled" changes don't conflict with your code (i.e. they are to files that are different to the one you've been editing). In that case, you don't need to do anything, just push your commit to the repo.

2. The "pulled" changes conflict with your file, but in a different place to the bit you've been editing. In this case git tries to "auto-merge" the two files together. The message will say "Auto-merging XXX". When this has updated, you should see edits that someone else has made to your code, along with your edits. You can then "push" your changes to the server.

3. The "pulled" changes conflict with your file in the same bit that you've been editing. In this case, git will add the conflict to the code, with a bunch of lines showing what your stuff looks like (labeled as "HEAD"), and the bits that someone else has been working on. You now need to talk to your friend who did the editing, find out what changes should be kept. Then you should add, commit and push to the server. You can also use tools to help you pick out which changes to keep, but I haven't figured this out yet. It should follow from:

git mergetool