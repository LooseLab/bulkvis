import textwrap

_help = "Output the citation for this tool and exit"
_cli = ()


def run(parser, args):
    cite = textwrap.fill(
        (
            "Alexander Payne, Nadine Holmes, Vardhman Rakyan, Matthew Loose, "
            "BulkVis: a graphical viewer for Oxford nanopore bulk FAST5 files, "
            "Bioinformatics, Volume 35, Issue 13, 1 July 2019, Pages 2193–2198"
        ),
        width=70,
        subsequent_indent=" " * 10,
    )
    url = "https://academic.oup.com/bioinformatics/article/35/13/2193/5193712"
    doi = "10.1093/bioinformatics/bty841"
    print("Thank you for using bulkvis!\n")
    print(f"Citation: {cite}")
    print(f"URL:      {url}")
    print(f"DOI:      {doi}")
