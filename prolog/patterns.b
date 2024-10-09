%search constraints

:- set(i,20).
:- set(clauselength,4).
:- set(search,ar).
:- set(pos_fraction,0.001).
:- set(max_features,4097).
:- set(noise,inf).
:- set(nodes,200000). % let the search run its course


% mode declarations (note that ORF and gene are interchangeable in this specific context)
:- modeh(*, gene(+orf)).
:- modeb(*, compound_name(+value, #compound)).
:- modeb(*, drug_target(+value, #action, -orf_through)).
:- modeb(*, condition(+cond, #cond_desc)).
:- modeb(*, phenotype(+orf, #state, -value, -cond)).
:- modeb(*, interacts_with_metabolite(+orf_through, #type, #reaction, #metabolite)).
:- modeb(*, perturbation_in_metabolism_of(+orf, #metabolite2, #type2, #reaction2)).

% determinations
:- determination(gene/1, compound_name/2).
:- determination(gene/1, drug_target/3).
:- determination(gene/1, condition/2).
:- determination(gene/1, phenotype/4).
:- determination(gene/1, interacts_with_metabolite/4).
:- determination(gene/1, perturbation_in_metabolism_of/4).

% loading in background knowledge
:- ['kb/phenotypes.dl'].
:- ['kb/drug_name.dl'].
:- ['kb/target.dl'].
:- ['kb/conditions.dl'].
:- ['kb/metabolites.dl'].

% Change the naming convention and order of the constants, as to make it more readable in natural language (for the meaning it is intended to give
% when having it directly connected through a positive example, or indirectly through a drug target).

perturbation_in_metabolism_of(Orf, Metabolite, Type, Reaction) :-
	interacts_with_metabolite(Orf, Reaction, Type, Metabolite).


% show examples as boolean vectors, included to provide an output of the "induce features" mode, as per the aleph manual. Allows for printing out and saving the features.
:- set(portray_examples,true).
aleph_portray(train_pos):-
        setting(train_pos,File),
        show_features(File,positive).
aleph_portray(train_neg):-
        setting(train_neg,File),
        show_features(File,negative).
aleph_portray(test_pos):-
        setting(test_pos,File),
        show_features(File,positive).

show_features(File,Class):-
        open(File,read,Stream),
        repeat,
        read(Stream,Example),
        (Example = end_of_file -> close(Stream);
                write_features(Example,Class),
                fail).
write_features(Example,_):-
        feature(_,(Example:- Body)),
        (Body -> write(1), write(' '); write(0), write(' ')),
        fail.
write_features(_,Class):-
	writeq(Class), nl.


:- ['saveOutput.pl']. % function to allow for saving the output
