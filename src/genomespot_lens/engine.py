"""Adapter to the pinned GenomeSPOT source; tested against its public API.

Per-protein metrics and DNA statistics are invariant to protein selection.
Localization and aggregate statistics are recalculated for EVERY subset.
The cache uses upstream private fields, so changing versions requires parity tests.
"""

from pathlib import Path
import subprocess

from . import UPSTREAM_COMMIT


class Engine:
    def __init__(self, contigs, proteins, models):
        import genome_spot
        from genome_spot.bioinformatics.genome import Genome
        from genome_spot.genome_spot import GenomeSPOT

        root = Path(genome_spot.__file__).resolve().parent.parent
        try:
            revision = subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"], text=True).strip()
            dirty = subprocess.check_output(["git", "-C", str(root), "status", "--porcelain", "--untracked-files=no"], text=True).strip()
        except (subprocess.CalledProcessError, FileNotFoundError) as exc:
            raise ValueError("Install GenomeSPOT editable from its pinned git checkout; see README") from exc
        if revision != UPSTREAM_COMMIT or dirty:
            raise ValueError("GenomeSPOT checkout differs from the tested, unmodified pinned revision")
        self.root = root
        self.models = str(Path(models).resolve())
        self.genome = Genome(str(contigs), str(proteins))
        self.metrics = self.genome.protein_data()
        self.dna = self.genome.compute_dna_statistics()
        self.predictor = GenomeSPOT()

    def features(self, ids):
        if not ids or not set(ids) <= self.metrics.keys():
            raise ValueError("Invalid selection for cached GenomeSPOT metrics")
        genome = self.genome
        genome._protein_data = {key: self.metrics[key] for key in sorted(ids)}
        genome._protein_localization = None
        locations = genome.assign_localization()
        features = {"all": dict(self.dna), "diff_extra_intra": {}}
        features["all"].update(genome.compute_protein_statistics())
        features["all"]["protein_coding_density"] = 3 * features["all"]["total_protein_length"] / self.dna["nt_length"]
        for category, upstream_label in [("extracellular_soluble", "extra_soluble"), ("intracellular_soluble", "intra_soluble"), ("membrane", "membrane")]:
            subset = {key for key, value in locations.items() if value == upstream_label}
            if subset:
                features[category] = genome.compute_protein_statistics(subset)
        # Missing compartments stay missing, yielding upstream missing-feature warnings.
        # Upstream raises KeyError for all-extracellular subsets; avoid that here.
        if "intracellular_soluble" in features:
            for key, value in features.get("extracellular_soluble", {}).items():
                features["diff_extra_intra"][key] = value - features["intracellular_soluble"].get(key, float("nan"))
        return features

    def predict(self, ids):
        features = self.features(ids)
        predictions = self.predictor.predict_from_genome(features, self.models)
        counts = {k: v.get("total_proteins", 0) for k, v in features.items() if k != "diff_extra_intra"}
        return predictions, counts
