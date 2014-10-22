------------------------------------------------------------------
Overview
------------------------------------------------------------------
This is the code that supports the SCFG grammar sampler described in 
"A Bayesian Model for Learning SCFGs with Discontiguous Rules" 
(http://aclweb.org/anthology/D/D12/D12-1021.pdf). As well it includes 
an extended version of the sampler that incorporates block sampling but
has strong dependencies to the cdec decoder (http://cdec-decoder.org/index.php?title=Main_Page)
to support the bilingual parsing. (Work is currently is submission.) 

------------------------------------------------------------------
Compiling the binary sampler 
------------------------------------------------------------------
Typing 
>> make binary
in the pyp-release directory should compile all the necessary dependencies. 

------------------------------------------------------------------
Compiling the block sampler 
------------------------------------------------------------------
For simplicity the current version of the sampler is released with the code for 
the block sampler commented out. If you would like to test the block sampler 
uncomment the following lines:
src/multiNTs.h lines 101
src/multiNTs.cpp lines 3, 516-519, 540, 552, 561  

Compile by typing
>> make block
in the pyp-release directory. Note that you need cdec installed and the compile
command must have access to all the header files and dependencies in the block_sampler.h
include list. See the Makefile for an example.  

------------------------------------------------------------------
Running the binary sampler 
------------------------------------------------------------------
After compiling run the sampler with a command such as the following:
>> ./multiNTs -s data/urdu.txt -t data/english.txt -m1 data/ur-en.t1.15 -m1i data/en-ur.t1.15 -m1a data/model1.gdfa -i 20 -s2d -1 -threads 2 -o outputfiles/toy

The parameters are
-s: source language text file with a sentence per line.
-t: target language text file with a sentence per line.
-m1: IBM Model1 translation probability file.
-m1i: IBM Model1 translation probability file.
-m1a: alignment initialisation.
-i: number of iterations to run for.
-s2d: acronym for "sentence to debug". Set to -1 to turn off sentence specific debugging.   
-threads: number of threads to use.
-o: is the output directory (./outputfiles/<unix-timestamp> by default). 

The full parameter descriptions are at the top of src/multiNTs.h.

------------------------------------------------------------------
Running the block sampler 
------------------------------------------------------------------
After compiling run the block sampler with a command such as 
>> ./multiNTs -s data/urdu.txt -t data/english.txt -m1 data/ur-en.t1.15 -m1i data/en-ur.t1.15 -m1a data/model1.gdfa -i 20 -s2d -1 -threads 2 -bs -bs-every 5

The parameters are
-bs: run the block sampler
-bs-every: how often to run the block sampling inbetween binary samping iterations

The rest of the parameters are the same as the binary sampler. 

------------------------------------------------------------------
Input files
------------------------------------------------------------------
data/en-ur.t1.15: IBM1 English->Urdu word translation probabilities
data/english.txt: English training corpus
data/model1.gdfa: IBM1 word alignments for PYP initialization
data/ur-en.t1.15: IBM1 Urdu->English word translation probabilities
data/urdu.txt: Urdu training corpus 

------------------------------------------------------------------
Output files
------------------------------------------------------------------
Files are saved with the following named conventions.

alignment.n: word alignments for iteration n. n = -1 indicates initialization files. 
rules.n: the PYP-SCFG grammar for iteration n. Again n = -1 indicates initialization. 
trees.state.gz: saves the current iteration's PYP-SCFG grammar state to restart sampling
                from the last iteration. Details are below. 

------------------------------------------------------------------
Restarting sampling
------------------------------------------------------------------
The file in the output directory trees.state.gz saves the current state of the PYP-SCFG grammar.
To restart sampling from the last iteration use a command such as the following:

>> ./multiNTs -s data/urdu.txt -t data/english.txt -m1 data/ur-en.t1.15 -m1i data/en-ur.t1.15 -m1a outputfiles/1375447423/trees.state.gz -restart -i 20 -s2d -1  -o outputfiles/toy/restart/

where
-restart: signals a restart 
-m1a: points to the initial alignments to use in the trees.state.gz file
-o: points to the output directory to save the restarted experiment files. 

The rest of the parameters are all the same as before. 
