"""Read one sequence/scaffold record and describe its supplied scaffold."""

from pathlib import Path


def scaffold_pairs(scaffold):
    stack, pairs = [], set()
    for position, symbol in enumerate(scaffold, 1):
        if symbol == "(":
            stack.append(position)
        elif symbol == ")":
            if not stack:
                raise ValueError("unbalanced scaffold")
            pairs.add((stack.pop(), position))
        elif symbol != ".":
            raise ValueError("scaffold must be pseudoknot-free dot-bracket")
    if stack:
        raise ValueError("unbalanced scaffold")
    return pairs


def paired_positions(pairs):
    return {position for pair in pairs for position in pair}


def read_record(path):
    lines = [line.strip() for line in Path(path).read_text().splitlines()
             if line.strip() and not line.lstrip().startswith("#")]
    if len(lines) != 3 or not lines[0].startswith(">"):
        raise ValueError("accuracy inputs require exactly one named record")
    name, sequence, scaffold = lines[0][1:], lines[1], lines[2]
    if not name or not sequence or set(sequence) - set("ACGU"):
        raise ValueError("a named ACGU sequence is required")
    if len(sequence) != len(scaffold):
        raise ValueError("sequence/scaffold length mismatch")
    pairs = scaffold_pairs(scaffold)
    if any(sequence[i - 1] + sequence[j - 1] not in
           {"AU", "UA", "CG", "GC", "GU", "UG"} for i, j in pairs):
        raise ValueError("scaffold contains a noncanonical pair")
    return name, sequence, scaffold


def describe_scaffold(sequence, scaffold):
    pairs = scaffold_pairs(scaffold)
    children = {None: []}
    stack = []
    maximum_depth = 0
    for position, symbol in enumerate(scaffold, 1):
        if symbol == "(":
            parent = stack[-1] if stack else None
            children[parent].append(position)
            children[position] = []
            stack.append(position)
            maximum_depth = max(maximum_depth, len(stack))
        elif symbol == ")":
            stack.pop()
    branching_depth = int(len(children[None]) >= 2)
    maximum_branching_ancestors = 0
    stack = []
    for position, symbol in enumerate(scaffold, 1):
        if symbol == "(":
            stack.append(position)
            branching_depth += int(len(children[position]) >= 2)
        elif symbol == ")":
            branching_depth -= int(len(children[stack.pop()]) >= 2)
        else:
            maximum_branching_ancestors = max(maximum_branching_ancestors, branching_depth)
    return {"length": len(sequence), "scaffold_pairs": len(pairs),
            "free_bases": len(sequence) - 2 * len(pairs),
            "paired_fraction": 2 * len(pairs) / len(sequence),
            "maximum_pair_depth": maximum_depth,
            "maximum_branching_ancestors_of_free_base": maximum_branching_ancestors,
            "branching_nodes_including_exterior":
            sum(len(nodes) >= 2 for nodes in children.values()),
            "gc_fraction": sum(base in "GC" for base in sequence) / len(sequence)}
