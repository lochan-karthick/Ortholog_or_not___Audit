# pylint: disable=W0231
from abc import ABC, abstractmethod
from collections import UserList
import dataclasses
from glob import glob
import os.path
import re
import yaml

import dask.dataframe as dd
from gffutils.exceptions import FeatureNotFoundError
from natsort import natsorted
import pysam
from tqdm.dask import TqdmCallback

from utils.gffutils import init_db
from .constants import BLAST_COLUMNS
from .utils import SpeciesIDMapping

WBPS_RELEASE = "WBPS19"


class SpeciesList(UserList):
    def __init__(self, initlist, **kwargs):
        wd_path = kwargs["wd_path"]
        load_blast = kwargs.get("load_blast", False)
        self.data = []
        for sp in initlist:
            sp.id = SpeciesIDMapping(wd_path, sp.prot_meta.filename_suffix)[sp.data_label]
            seen_pairs = []
            blast_paths = []
            for bp in natsorted(glob(os.path.join(wd_path, "Blast*"))):
                pair = sorted(map(int, re.findall(r'\d+', os.path.basename(bp))))
                if sp.id in pair and pair not in seen_pairs and len(set(pair)) == 2:
                    if load_blast:
                        print(f"loading {bp}...")
                    blast_paths.append(bp)
                    seen_pairs.append(pair)
            if load_blast:
                blastout = dd.read_csv(blast_paths, names=BLAST_COLUMNS, delimiter="\t", dtype={'bitscore': 'float64'})
                sp.slice_blast_data(blastout)
            self.data.append(sp)

    def get_species_ids_for_clade(self, clade):
        species_ids = set()
        for sp in self:
            if sp.clade == clade:
                species_ids.add(str(sp.id))
        return tuple(species_ids)

    def get_species_with_data_label(self, data_label):
        for sp in self:
            if sp.data_label == data_label:
                return sp


@dataclasses.dataclass
class ProteinMeta:
    file_path: str
    filename_suffix: str
    label: str


class AltSourceMixin:
    release = ""

    def __init__(self, *args, **kwargs):
        self.data_dir = kwargs.pop("data_dir", self.data_dir)
        super().__init__(*args, **kwargs)

    @property
    def db_path(self):
        return os.path.join(self.db_dir, f"{self.abbr}{self.name}.db")

    @property
    def gff_path(self):
        return os.path.join(self.data_dir, ".".join((self.data_label, "gff3")))


class Species(ABC):
    id: int
    tmp_bed_path: str
    data_dir = os.path.join("data", "from_WBPS", "")
    db_dir = "db"
    release = WBPS_RELEASE
    gff_id_prefix = "transcript:"
    hog_id_prefix = "transcript_"
    sequence_id_prefix = "transcript:"
    
    def get_gff_transcript_id(self, tid):
        if self.gff_id_prefix and tid.startswith(self.gff_id_prefix):
            return tid

        return self.gff_id_prefix + tid
    
    def get_internal_transcript_id(self, tid):
        if self.gff_id_prefix and tid.startswith(self.gff_id_prefix):
            return tid[len(self.gff_id_prefix):]

        return tid
    
    def get_hog_transcript_id(self, tid):
        if self.hog_id_prefix and tid.startswith(self.hog_id_prefix):
            return tid[len(self.hog_id_prefix):]

        return tid
    
    def get_sequence_transcript_id(self, tid):
        if self.sequence_id_prefix and tid.startswith(self.sequence_id_prefix):
            return tid[len(self.sequence_id_prefix):]

        return tid
    
    def get_gene_ids(self, tid):
        gff_tid = self.get_gff_transcript_id(tid)

        gene_ids = {
            parent.id
            for parent in self.db.parents(
                gff_tid,
                featuretype="gene"
            )
        }

        if not gene_ids:
            raise ValueError(
                f"No gene parent found for transcript {gff_tid}"
            )

        return gene_ids
    
    
    default_transcript_selection_method = "_get_longest_transcript"
    

    def __init__(self, name, acc="", prot_filename_suffix=".protein.fa", data_label=None, skip_plot=False):
        self.name = name
        self.prefix = self.abbr + self.name.lower()[:3]
        self.data_label = f"{self.genus}_{self.name}.{acc}.{self.release}" if not data_label else data_label
        self.db = init_db(self.gff_path, self.db_path)
        self.blast_slice = None
        self.prot_meta = ProteinMeta(
            file_path=os.path.join(self.data_dir, self.data_label + prot_filename_suffix),
            filename_suffix=prot_filename_suffix,
            label=self.data_label + ".".join(prot_filename_suffix.split(".")[:-1])
        )
        self.all_transcript_ids = [line.strip(">").strip("transcript:").split(" ")[0] for line in open(self.prot_meta.file_path) if line.startswith(">")]
        self.skip_plot = skip_plot

    @property
    def gff_path(self):
        return os.path.join(self.data_dir, ".".join((self.data_label, "annotations", "gff3")))

    @property
    def db_path(self):
        return os.path.join(self.db_dir, f"{self.prefix}_wbps.db")

    @property
    @abstractmethod
    def abbr(self):
        raise NotImplementedError()

    @property
    @abstractmethod
    def genus(self):
        raise NotImplementedError()

    def slice_blast_data(self, blastout):
        left = blastout[(blastout["qseqid"].str.startswith(str(self.id))) & ~(blastout["sseqid"].str.startswith(str(self.id)))]
        left = left.rename(columns={"qseqid": "transcript_id", "sseqid": "other_transcript_id"})
        right = blastout[(blastout["sseqid"].str.startswith(str(self.id))) & ~(blastout["qseqid"].str.startswith(str(self.id)))]
        right = right.rename(columns={"qseqid": "other_transcript_id", "sseqid": "transcript_id"})
        df = dd.concat([left, right], axis=0)
        self.blast_slice = df

    def select_transcript(self, prot_ids:set, other_prot_ids:set, seq_id_map, load_blast=False):
        if load_blast and self.blast_slice is not None and other_prot_ids:
            if hasattr(self.blast_slice, "compute"):
                with TqdmCallback(desc=f"compute blast slice for species {self.id}"):
                    self.blast_slice = self.blast_slice.compute()
            return self._get_transcript_with_best_blast(prot_ids, other_prot_ids, seq_id_map)
        return getattr(self, self.default_transcript_selection_method)(prot_ids)

    def _get_transcript_with_most_exons(self, transcript_ids):
        exon_counts = {}
        
        for tid in transcript_ids:
            gff_tid = self.get_gff_transcript_id(tid)
            cds = list( self.db.children( gff_tid, featuretype="CDS"))
            exon_counts[len(cds)] = { self.db[gff_tid]: cds }
        for tran, exons in exon_counts[max(exon_counts)].items():
            return tran, sorted(exons, key=lambda x: x.start, reverse=tran.strand == "-")

    def _get_longest_transcript(self, transcript_ids):
        prot_lengths = {}
        for tid in transcript_ids:
            gff_tid = self.get_gff_transcript_id(tid)
            cds = list(self.db.children(gff_tid,featuretype="CDS"))
            prot_lengths[self.get_amino_acid_count(cds)] = {self.db[gff_tid]: cds}
        for tran, exons in prot_lengths[max(prot_lengths)].items():
            return tran, sorted( exons, key=lambda x: x.start, reverse=tran.strand == "-")    

    def _get_transcript_with_best_blast(self, transcript_ids, other_transcript_ids, seq_id_map):
        filt = self.blast_slice[(self.blast_slice["transcript_id"].isin(map(seq_id_map.get, transcript_ids))) &
                                (self.blast_slice["other_transcript_id"].isin(map(seq_id_map.get, other_transcript_ids)))]
        tid = seq_id_map[filt.groupby("transcript_id").median(["pident", "bitscore"]).sort_values(["pident", "bitscore"], ascending=False).head(1).index[0]]
        try:
            tran = self.db[tid]
        except FeatureNotFoundError:
            tran = self.db[self.get_gff_transcript_id(tid)]
        exons = list(self.db.children(self.get_gff_transcript_id(tid) ,featuretype="CDS"))
        return tran, sorted(exons, key=lambda x: x.start, reverse=tran.strand=="-")

    def get_protein_sequence(self, tid):
        # Slow method to obtain amino acids using pysam...
        seq = pysam.faidx(self.prot_meta.file_path, tid)
        return "".join(seq.split("\n")[1:])

    @staticmethod
    def get_amino_acid_count(cds_exons):
        # ...or faster method by inferring from CDS lengths
        return int(sum(len(cds) for cds in cds_exons) / 3)
    
    
#Describes species whose metadata is supplied through the species config file.
class ConfiguredSpecies(Species):
    
    def __init__(
        self,
        name,
        genus,
        abbr,
        clade,
        data_dir,
        data_label,
        gff_filename,
        prot_filename_suffix=".fa",
        gff_id_prefix="",
        hog_id_prefix="",
        sequence_id_prefix="",
        skip_plot=False
    ):
        self._genus = genus
        self._abbr = abbr
        self.clade = int(clade)
        self.data_dir = data_dir
        self.gff_filename = gff_filename
        self.gff_id_prefix = gff_id_prefix
        self.hog_id_prefix = hog_id_prefix
        self.sequence_id_prefix = sequence_id_prefix
        
        super().__init__(
            name,
            data_label=data_label,
            prot_filename_suffix=prot_filename_suffix,
            skip_plot=skip_plot
        )
    @property
    def genus(self):
        return self._genus

    @property
    def abbr(self):
        return self._abbr

    @property
    def gff_path(self):
        return os.path.join(
            self.data_dir,
            self.gff_filename
    )

    @property
    def db_path(self):
        return os.path.join(
            self.db_dir,
            f"{self.abbr}{self.name}.db"
        )


class Schistosoma(Species):
    abbr = "S"
    genus = "schistosoma"


class HaematobiumClade(Schistosoma):
    clade = 1


class HaematobiumCladeFromTool(AltSourceMixin, HaematobiumClade):
    pass


class MansoniClade(Schistosoma):
    clade = 2


class MansoniCladeFromTool(AltSourceMixin, MansoniClade):
    pass


class JaponicumClade(Schistosoma):
    clade = 3


class JaponicumCladeFromTool(AltSourceMixin, JaponicumClade):
    pass

def load_species_config(config_path, data_dir):
    species = []

    with open(config_path, "r") as config_file:
        config = yaml.safe_load(config_file)

    if not isinstance(config, dict) or "species" not in config:
        raise ValueError(
            "Species config must contain a top-level 'species' list"
        )

    if not isinstance(config["species"], list):
        raise ValueError(
            "'species' must be a list"
        )

    required_fields = {
        "name",
        "genus",
        "abbr",
        "clade",
        "data_label",
        "prot_filename_suffix",
        "files",
        "id_prefixes",
        "skip_plot"
    }

    required_prefixes = {
        "hog",
        "sequence",
        "gff"
    }

    for entry in config["species"]:
        missing = required_fields.difference(entry)

        if missing:
            raise ValueError(
                "Species config entry is missing required fields: "
                + ", ".join(sorted(missing))
            )
            
        files = entry["files"]

        if not isinstance(files, dict):
            raise ValueError(
                f"files must be a mapping for species {entry['name']}"
            )

        if "gff" not in files:
            raise ValueError(
                f"Species {entry['name']} is missing GFF file information"
            )

        prefixes = entry["id_prefixes"]

        if not isinstance(prefixes, dict):
            raise ValueError(
                f"id_prefixes must be a mapping for species {entry['name']}"
            )

        missing_prefixes = required_prefixes.difference(prefixes)

        if missing_prefixes:
            raise ValueError(
                f"Species {entry['name']} is missing ID prefixes: "
                + ", ".join(sorted(missing_prefixes))
            )

        try:
            clade = int(entry["clade"])
        except (TypeError, ValueError):
            raise ValueError(
                f"clade must be an integer for species {entry['name']}"
            )

        skip_plot = entry["skip_plot"]

        if not isinstance(skip_plot, bool):
            raise ValueError(
                f"skip_plot must be true or false for species {entry['name']}"
            )

        species.append(
            ConfiguredSpecies(
                name=entry["name"],
                genus=entry["genus"],
                abbr=entry["abbr"],
                clade=clade,
                data_dir=data_dir,
                data_label=entry["data_label"],
                gff_filename=files["gff"],
                prot_filename_suffix=entry["prot_filename_suffix"],
                hog_id_prefix=prefixes["hog"],
                sequence_id_prefix=prefixes["sequence"],
                gff_id_prefix=prefixes["gff"],
                skip_plot=skip_plot
            )
        )

    return species


class IndicumClade(Schistosoma):
    clade = 4


class NewSpeciesClade(Schistosoma):
    clade = 5


class Haemonchus(Species):
    abbr = "H"
    genus = "haemonchus"
    clade = 0


class HaemonchusFromTool(AltSourceMixin, Haemonchus):
    pass


class Pristionchus(Species):
    abbr = "P"
    genus = "pristionchus"
    clade = 0


class PristionchusFromTool(AltSourceMixin, Pristionchus):
    pass
