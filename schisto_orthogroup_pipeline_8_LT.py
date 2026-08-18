#! /usr/bin/env python3

from orthologue_analysis import (
    parse_args,
    main
)
from orthologue_analysis.species import (
    SpeciesList,
    load_species_config,
)

if __name__ == "__main__":
    args = parse_args()
    args.prefix_cut = "transcript_"

    # Load the species defined by the user in the configuration file
    species = load_species_config(
        args.species_config,
        args.species_data_dir
    )

    SPECIES_LIST = SpeciesList(
        species,
        **vars(args)
    )

    main(args, SPECIES_LIST)
