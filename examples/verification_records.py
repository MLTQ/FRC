"""Readable object records with compact numeric arrays, without data changes."""
import json


def _array_tree(value):
    if isinstance(value, dict):
        return False
    return all(_array_tree(item) for item in value) if isinstance(value, list) else True


def _render(value, depth=0):
    indent = "  "*depth
    if isinstance(value, dict):
        if not value:
            return "{}"
        return "{\n"+",\n".join("  "+indent+json.dumps(key)+": "+_render(item, depth+1)
                                  for key, item in value.items())+"\n"+indent+"}"
    if isinstance(value, list) and not _array_tree(value):
        return "[\n"+",\n".join("  "+indent+_render(item, depth+1) for item in value)+"\n"+indent+"]"
    return json.dumps(value, separators=(",", ":"), allow_nan=False)


def write_record(path, value):
    text = _render(value)+"\n"
    if json.loads(text) != value:
        raise ValueError("Scientific record serialization changed its data")
    path.write_text(text)
