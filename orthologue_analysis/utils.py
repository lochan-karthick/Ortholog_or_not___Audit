import os.path

from utils.generic import get_project_root


class SequenceIDMapping:
    def __init__(self, wd_path, species_list):
        self.map = {}
        self.inv_map = {}

        #Match each OrthoFinder species ID to its Species object
        species_by_id = {
            str(sp.id): sp
            for sp in species_list
        }

        with open(os.path.join(wd_path, "SequenceIDs.txt"), "r") as f:
            for l in f:
                sid, info = l.strip().split(": ")
                tid = info.split(" ")[0]

                #OrthoFinder IDs are formatted as speciesID_sequenceID, for example 6_9762.
                species_id = sid.split("_", 1)[0]

                #Ignore species that are not part of this analysis
                sp = species_by_id.get(species_id)

                if sp is None:
                    continue

                #Convert the SequenceIDs.txt transcript ID into the pipeline's standard internal transcript ID.
                
                tid = sp.get_sequence_transcript_id(tid)

                self.map[tid] = sid
                self.inv_map[sid] = tid

    def __getitem__(self, item):
        item = str(item)

        if item in self.inv_map:
            return self.inv_map[item]

        return self.map[item]

    def get(self, item):
        return self[item]


class SpeciesIDMapping:
    def __init__(self, wd_path):
        self.map = {}

        with open(os.path.join(wd_path, "SpeciesIDs.txt"), "r") as f:
            for l in f:
                sid, prot_path = l.strip().split(": ")

                protein_filename = os.path.basename(prot_path)

                self.map[protein_filename] = int(sid)

    def __getitem__(self, item):
        return self.map[os.path.basename(item)]


def orthofinder_paths(label, subdir = "Phylogenetic_Hierarchical_Orthogroups"):
    orthofinder_dir = get_project_root() / "data" / "from_MARS" / "OrthoFinder"
    paths = {}
    paths["wd"] = orthofinder_dir / "WorkingDirectory" / label
    if subdir == "Phylogenetic_Hierarchical_Orthogroups":
        paths["orthogroups"] = orthofinder_dir / subdir / label / "N0.tsv"
    elif subdir == "Orthogroups":
        paths["orthogroups"] = orthofinder_dir / subdir / label / "Orthogroups.tsv"
        paths["orthogroups_unassigned_genes"] = orthofinder_dir / subdir / label / "Orthogroups_UnassignedGenes.tsv"
    return paths
