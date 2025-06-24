cell(A) :-
    exhibits_phenotype(A,increased, "resistance to chemicals",B,C),
    compound_name(B,fenpropimorph),
    condition(C,normal).

exhibits_phenotype(A,D,"resistance to chemicals",B,C):-
    triple(o_001,isa,resistanceToChemicals).