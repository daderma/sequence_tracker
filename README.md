## Getting the latest version of the program ##

Go to the "Releases" tab and download latest release archive. It is a 64-bit Windows executable that will run on any recent Windows version. 

In addition the program can be built from source for Mac and Linux.


## How to extract data from immunoSEQ Analyzer ##

1) Include the desired samples in the project.

2) Create a new analysis containing all samples.

3) Open the analysis in "Advanced Query" and run this query: "**select * from sequences**"

4) Save the output to file (must be .tsv extension).


## Performing iterative analysis (tracking) ##

To track sequences across samples (could be different tissue types or over time) you must tag the sample data.

Each child sample that you want to compare to a parent sample must contain a tag indicating the name of the parent. A child sample can contain multiple parents, if you want to branch the analysis.

1) Choose the appropriate project in immunoSEQ Analyzer and click "Analyze".

2) In project overview, start by selecting/copying (ctrl-C) the sample name of the parent.

3) Select the child sample and choose "Edit tags".

4) In the right hand panel, choose "Create a project-specific tag".

5) In "Tag Group" select "Create new tag group" (you only need to do this the first time). Name of the group must be "**ParentSample**".

6) Next, create the new tag with "Tag Group" set to "ParentSample" and paste the name of the parent sample into the "Tag" value.

7) Extract the data as explained above.
