# A Bunch of thoughts

---
## sAMPC sucks, like legit

sAMPC kinda sucks, and lets be honest the code is un-maintainable for someone other than me at the moment.
Nevertheless, you may have noticed that some of the code on the development branch has `Type annotations`.
I've started implementing this to make it easier to read and will be implementing it to all the code before
merging it with the main branch.

You may have also noticed that the code is broken down into multiple files on the development branch, again this
is an attempt by me to make the code easier to work with.

### A Few horrors
1. Apparently in the feature extraction functions I have merged columns from DataFrames and just straight up dropped indexes. 
Yeah....**DROPPED THE INDEXES**, I am working on a fix for that right now, but apparently it works correctly. 
Please do not trust that no matter what.

2. The base models used for the stack are tested based on a subset of the training dataset. This would be fine normally, due to 
the limited amount of experimentally verified true positives, there is a significant overlap between the two datasets.
I don't really know how to fix this honestly because it is not in my hands, I could do a K-fold cross validation but since 
I am training literally thousands of models I have no idea how it will affect my PC. One other option I have is to use predicted 
AMPs from the AMPSphere Database, which could be useful, but then again, I got no idea how that would affect the models and what biases 
these models would have.

3. There are absolutely no sanity checks in place, which makes it really annoying to work with if you don't know what to look for. 
Few DataFrames are saved as Excel files and have hard coded columns, so if someone messes up the column names there is no way to figure out 
what went wrong.

4. There are no DocStrings.

Here is a to-do list on the stuff I am planning to implement
- [x] Fix the logic bug
- [] Add DocStrings for all the functions to make it maintainable 
- [] Fix a few horrors
- [] Check if everything I did works with `pandas 3.0`
- [] Figure out a way to do the base model testing without relying on the training dataset
- [x] Fix the scaling issue and figure out a way to be able to use the correct scaler for the correct Neural Network
- [] Improve on the fleiss_kappa(), currently it only takes binary inputs, 1 or 0


## Stuff I'd like to do to sAMPC in the future
1. I've implemented my own rendition of Fliess' Kappa in a function called fleiss_kappa() in `Base_model_selection.py`.
Currently, this is limited to binary inputs, 1 and 0 which is kinda annoying because then I needed to write my own logic to 
categorize my continuous predictions into True or False based on if they're above 0.5 or below 0.5. Also, if I need to re-use this function 
for something else or someone decides to play with the model agreement/disagreement by using different categories, they'd need to write their 
own logic since the binary logic is hard coded into the function. To fix this, I'd like to implement some way to categorize the data based on 
number of categories defined by a user or something that can find the best number of categories for that data. SKlearn's Random Forest classifier 
uses something similar to find the best way the data should be split for its trees if I am not mistaken. So maybe I'd reverse engineer that, idk its 
just a thought.

2. I want to make sAMPC installable using pip, I am new to this so I got no idea how to do this yet. I'd also like to make the standalone tool that can 
run like any bioinformatics tool with its own conda enviornment, but with the diabolical amount of user inputs each function has, it would be a challenge.

3. I really want to try using DNA Sequences, Protein Sequences and embeddings from those fancy models by Google and META. I have no idea how to do this yet 
or what the plan is but since sAMPC is also a learning project, why not try it out?

## The bug that keeps getting worse the more you think about it
So, if you've read the changelog on the development branch then you've seen that there was a logical bug in the stack_helper().
At first glance this looks like a simple innocent little indentation mistake, which it is -- an extra indent that puts the averaging step 
inside the for loop. So, the function takes an average every time a column (here the predictions of a model on the test set) is added to the DataFrame.
It is still not a huge deal, because if you check the math, it is essentially adding a weight, where low scoring proteins get a lower average and high scoring 
proteins get a higher average. 
If you look closely, you'll also notice that this "weight" is mostly dependent on the first column's (first model's) prediction. If the model says the protein is 1 
then the first average starts with 1, if the model says the protein is 0, then the first average starts with 0. Okay....a little shady right? Now any sane person would
ask, "who selects which model comes first?"....and that my friend is where it starts getting worse, I haven't really implemented anything to select a model first, so the model
is loaded based on whatever way the user's device saves the file (by date, by name, ascending, descending, etc) which is essentially random. So, the results I got on my PC
may not reflect what a user will see on their device. Yaay production code (┛◉Д◉)┛彡┻━┻.

---

I know no one is really reading this but hey, it's nice to write these thoughts down for the future.
