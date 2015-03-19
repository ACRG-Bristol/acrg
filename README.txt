Shared Python code for the Atmospheric Chemistry Research Group (University of Bristol)

*First Installation*
For ACRG group members, check out the code using your (personal) Google account. Create a directory, and then, from that directory, type:

git clone https://mrghg@bitbucket.org/mrghg/acrg.git

*Adding a file to the repo*
From your directory, first add a file to the SVN library:

git add YOUR_FILE_NAME

Then "commit" all the changes you've made:

git commit -m "ADD_SOME_DESCRIPTIVE_COMMENTS_HERE_IN_QUOTES"

Note that it won't work if you don't add the -m with some comments.

Then "push" it to the central repository

git push

*Update regularly*
To make sure you're up-to-date with everyone else's changes, regularly type (or set a cron job), from inside the folder:

git checkout

*Changing someone else's code*
First, make sure you have the most up-to-date version of the repo (see above). Then, edit the code on your machine. Then, when you're done, just commit the changes back to the repo:

git commit -m "ADD_SOME_DESCRIPTIVE_COMMENTS_HERE_IN_QUOTES"

git push

*Adding an ssh key*

Git will ask you for a password every time you want to connect using https. You can get around this by adding your RSA key to your account. Instructions are here: https://confluence.atlassian.com/display/BITBUCKET/Set+up+SSH+for+Git. You can probably start at step 6.

