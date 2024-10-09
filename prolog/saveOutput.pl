%saveQueries(FileName) :-
%    protocol(FileName).
%
%stopQueriesSaving :-
%    noprotocol.

saveQueries(FileName) :-
    term_string(FileName, FileNameString),
    protocol(FileNameString).

stopQueriesSaving :-
    noprotocol.