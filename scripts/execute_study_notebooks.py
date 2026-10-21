"""Execute the study notebooks with this Python interpreter and validate outputs."""
import json
from pathlib import Path
import sys

import nbformat
from nbclient import NotebookClient
from jupyter_client import KernelManager
from jupyter_client.kernelspec import KernelSpecManager

root = Path(__file__).resolve().parents[1]
kernel_root = root / "runs/amp-study/kernels"
kernel = kernel_root / "genomespot-study"
kernel.mkdir(parents=True, exist_ok=True)
(kernel / "kernel.json").write_text(json.dumps({"argv": [sys.executable, "-m", "ipykernel_launcher", "-f", "{connection_file}"], "display_name": "GenomeSPOT Study", "language": "python"}))
paths = [Path(p) for p in sys.argv[1:]] or sorted((root / "notebooks").glob("0*.ipynb"))
for path in paths:
    print("Executing", path.name, flush=True)
    notebook = nbformat.read(path, as_version=4)
    manager = KernelManager(kernel_name="genomespot-study", kernel_spec_manager=KernelSpecManager(kernel_dirs=[str(kernel_root)]))
    client = NotebookClient(notebook, km=manager, timeout=600, resources={"metadata": {"path": str(root)}})
    try:
        client.execute()
    finally:
        if manager.has_kernel:
            manager.shutdown_kernel(now=True)
    nbformat.validate(notebook)
    errors = [output for cell in notebook.cells if cell.cell_type == "code" for output in cell.outputs if output.output_type == "error"]
    if errors:
        raise RuntimeError(f"Notebook contains execution errors: {path}")
    # Generic Python kernel metadata keeps the delivered notebooks portable.
    notebook.metadata.kernelspec.name = "python3"
    nbformat.write(notebook, path)
    print("Passed", path.name, flush=True)
