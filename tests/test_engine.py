import gzip
import os
from pathlib import Path

import numpy as np
import pytest

from genomespot_lens.core import read_fasta


@pytest.mark.integration
def test_cached_adapter_matches_upstream_full_and_subset(tmp_path):
    location = os.environ.get("GENOMESPOT_ROOT")
    if not location:
        pytest.skip("Set GENOMESPOT_ROOT to the pinned upstream checkout")
    from genome_spot.bioinformatics.genome import measure_genome_features
    from genome_spot.genome_spot import GenomeSPOT
    from genomespot_lens.engine import Engine

    root = Path(location)
    contigs = root / "tests/test_data/GCA_000172155.1_ASM17215v1_genomic.fna.gz"
    proteins = root / "tests/test_data/GCA_000172155.1_ASM17215v1_protein.faa.gz"
    records = read_fasta(proteins)
    engine = Engine(contigs, proteins, root / "models")
    for ids in [sorted(records), sorted(records)[::2]]:
        path = tmp_path / "subset.faa"
        path.write_text("".join(f">{k}\n{records[k]}\n" for k in ids))
        expected = measure_genome_features(str(path), str(contigs))
        actual = engine.features(ids)
        assert actual.keys() == expected.keys()
        for group, features in expected.items():
            assert actual[group].keys() == features.keys()
            for key, value in features.items():
                np.testing.assert_allclose(actual[group][key], value, rtol=1e-10, atol=1e-10, equal_nan=True)
        predictions, _ = engine.predict(ids)
        reference = GenomeSPOT().predict_from_genome(expected, str(root / "models"))
        for target, value in reference.items():
            for field, expected_value in value.items():
                if isinstance(expected_value, (float, np.floating)):
                    assert predictions[target][field] == pytest.approx(expected_value)
                else:
                    assert predictions[target][field] == expected_value
