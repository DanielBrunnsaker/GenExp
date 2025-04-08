- [ ] Rewrite Fuseki database script to write to `Dataset` rather than `Graph`.
- [ ] Implement logic for all phenotypes in Daniel's logic formulae
    - [ ] Resistance to chemicals

## Future

- [ ] If there are overlapping triples between hypotheses, write extra quads to the graph to include the triples in both (or indeed all) contexts. Currently, no new triples are written if they are all in the graph already. This can be problematic for reproducing versions of the dataset with different subsets of hypotheses. We need each hypothesis to be self-contained.
    - [ ] Will introduce duplication, maybe there is a neater way of doing this. Probably not though, without using something like RDF-star.
    - [ ] Need good queries for returing dataset versions given a subset of the hypotheses.