# Query, Fetch, Clean subcommand refinment/revamp

In there current form, Query, Fetch, Clean, Extract do not follow the conventions of the other tools of the phorge suite. This spec sheets hopes to corrct the flags, logic, and use of these subcommands.

## Overall notes

These subcommands in there current form share one thing. That is the flags do not follow the conventions set by the other subcommands in the tool kit, namely the user of flags and of the placement of inputs. This should be corrected and the commands should also be enhanced with features that will improve usability and functionality.

One area of improvement shared across all of these subcommands is that most use a form of json as input or output. I think it would be in the best intrest of the tool and of the user to transfer the files that are in json into TSV. I think this would allow for easier parsing of those files by a user. For example, a user can take a query file in tsv format, and using simple bash commands filter it to there own specification before passing it onto the fetch subcommand. 

second, that all log content should be written to stderr. some of the subcommands in this toolkit have required log output, such as concat. However, with these tool having the information about the progress of the subcommand will allow both the user to know what is going on while the subcommand runs and if they wish to record it they can simply capture it with ">>"  and modify it to whatever needs they have. This will be a more pragmatic approach that enables users more control of over information retention.

## query

```
$ phorge query -e you@email.com -t "COX1","12S" -i Pantherinae -o Hyaenidae,Felinae -q felidae_query.tsv
```
-e, --email: email address for NCBI API access
-i, --ingroup: ingroup taxa of intrest either a taxID for NCBI or Name (will need confirm prompt when using non taxaid query)
-o, --outgroup: outgroup taxa(s) either a taxID or Name
-q, --query: output query file in TSV format (same columns as json currently).
-t, --term: searh terms for loci or genetic data.
--api-key: API key for NCBI API access

Aside from the json output format becoming TSV as mentioned above, the use of A. search terms and B. common taxonomic names would help improve the use of this subcommand. Allowing for search terms to be used will enable the user to refine what they are looking for and not just pull everything as for some species that would be an impractical amount of data for most researchers questions. The use of being able to search with Common taxonomic names will improve the user interface (easier to remember a name, not a number) and readablility for reviewrs and collaboraters. With this feature will have to come a validation  of some kind such as "Found TAXON (taxID). is this correct [y/N]?". I could see potentially making the output go to stdout, but I fear that could lead to more confusion than is necessary. however, I am open to suggestions.

## fetch 

```
$ phorge fetch -e you@email.com  --min-length 100 --max-length 1000 -q felidae_query.tsv -o data/raw
```

-e, --email: email address for NCBI API access
--max-length: maximum length of sequences to fetch
--min-length: minimum length of sequences to fetch
-q, --query: input query file in TSV format (same columns as json currently)
-o, --output: output directory for fetched sequences
--api-key: API key for NCBI API access

improvements to the fetch subcommand come mostly form the TSV format change. In that, a user to able to parse and fileter the query TSV according to there needs, and the fetch subcomand is able to pull the data, as that as long as the accesion numbers are preserverd, data can be pulled. Further, if a user curates there own list of accesions but other means, they may simply pass a list of accesion numbers to fetch, increasing the interoperability of the tool with other common command line tools for bioinformatics. 

## clean

```
$ phorge clean -q felidae_query.tsv data/raw/*.fasta -o data/clean
```
-q, --query: input query file in TSV format as it acts as join table (same columns as json currently)
-o, --output: output directory for cleaned sequences
--prefer: way to allow preference over labs own voucher
-e, --extension: file extension for cleaned sequences "default: _std"

Again, mostly improvement comes from the TSV format change. However, includeing an extesion of some kind "_std","cln" would be helpful for tracking the progress of the data. 

## extract

```
$ phorge extract -r config/refrence.fna -m 0.7 -s 5.7 -c .5 data/clean/*.fna -o data/extracted/
```
-r, --reference: reference for of what you want extracted
-m, --min-identity: minimum identity threshold to be hit
-s, --sensitivity: sensitivity. Higher = more hits
-c, --coverage: minimum coverage threshold to be hit
--keep-intermediate: keep intermediate files from mmeqs2
--flank: number of bases to grab from both sites of the hit [default: 0]

Leverageing mmseqs is a great assest for this too. However other software that does something similer does exists and if there is a more pragmatic option, I would be willing to explore. A flag that should be added is -c (--coverage) which is not in the currect spec and is a attribute that is nessesary for making just judgement is extraction was successful and captured homologius sequences correctly.

# Conclusion

These subcommands are great. They serve very well. The were however a port from another tool and the flag conventions followed and were not updated to better integrate with the rest of the tool. These notes, along with discussion will bring them up to speed and enable there continued use in bioinformatics workflows.
