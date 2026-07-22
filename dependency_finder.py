from pathlib import Path
import ast
import json


OUTPUT_FILE = Path("dependency_audit.txt")


def find_python_imports():
    python_files = (
        sorted(Path(".").glob("schisto_orthogroup_pipeline*.py"))
        + sorted(Path("orthologue_analysis").glob("*.py"))
        + sorted(Path("utils").glob("*.py"))
    )

    output = []

    output.append("PYTHON IMPORTS")
    output.append("=" * 60)

    for path in python_files:
        output.append(f"\n[{path}]")

        try:
            tree = ast.parse(path.read_text(encoding="utf-8"))
        except Exception as error:
            output.append(f"Could not parse: {error}")
            continue

        imports = []

        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imports.extend(alias.name for alias in node.names)

            elif isinstance(node, ast.ImportFrom):
                module = "." * node.level + (node.module or "")
                names = ", ".join(alias.name for alias in node.names)
                imports.append(f"{module}: {names}")

        if imports:
            for imported_module in sorted(set(imports)):
                output.append(f"  {imported_module}")
        else:
            output.append("  No import statements found")

    return output


def find_notebook_imports():
    output = []

    output.append("\n\nNOTEBOOK IMPORTS")
    output.append("=" * 60)

    notebooks = sorted(
        Path(".").glob("schistosome_orthologue_analysis*.ipynb")
    )

    for path in notebooks:
        output.append(f"\n[{path}]")

        try:
            notebook = json.loads(path.read_text(encoding="utf-8"))
        except Exception as error:
            output.append(f"Could not read notebook: {error}")
            continue

        imports = []

        for cell in notebook.get("cells", []):
            if cell.get("cell_type") != "code":
                continue

            source = "".join(cell.get("source", []))

            for line in source.splitlines():
                line = line.strip()

                if line.startswith(
                    ("import ", "from ", "%run ", "!python")
                ):
                    imports.append(line)

        if imports:
            for import_line in sorted(set(imports)):
                output.append(f"  {import_line}")
        else:
            output.append("  No import statements found")

    return output


def main():
    output = find_python_imports() + find_notebook_imports()
    report = "\n".join(output)

    print(report)

    OUTPUT_FILE.write_text(report + "\n", encoding="utf-8")

    print(f"\nDependency report saved to: {OUTPUT_FILE.resolve()}")


if __name__ == "__main__":
    main()
